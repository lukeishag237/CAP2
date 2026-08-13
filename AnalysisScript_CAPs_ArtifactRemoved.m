% AnalysisScript_CAPs_ArtifactRemoved.m
%
% Preprocesses CAP data with artifact removal applied BEFORE epoch slicing.
% This script:
% 1. Loads raw data from LabChart files
% 2. Applies 60Hz filtering to the FULL raw sweep using filter_periodic_artifact
% 3. Slices epochs from filtered data
% 4. Applies artifact removal using template subtraction (from lowest amp)
% 5. Averages cleaned epochs
% 6. Saves to preprocessed_CAP_data_artifact_removed folder
%
% Uses the same approach as build_artifact_removed_df_SN.m for artifact removal.

function AnalysisScript_CAPs_ArtifactRemoved(expdate, file_prefix_adicht, file_prefix_outputs,...
                    CAP_avg_type, subtractnoise, zerobaseline,...
                    digitalLPF, digitalLPF_fc_bandpass, digitaloffsetremoval, ch_num_Vinput,...
                    ch_num_neural, stim_rec_dist, gain_value, fs, dt,...
                    labchart_unit_conv, stim_PW, stim_freq, V2mA, Ntrig_expected, phase_type,...
                    t_start_sig, CV_bounds, CV_bounds_slowest, y_min, y_max, y_max_f,...
                    amp_values_nominal, amp_values_actual, trial_values)

scrnsz      = get(0,'screensize');
indiv_plots = 0;

%% Add path to custom denoising scripts
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'external_dependencies', 'custom_denoising_scripts'));

%% Define data targets
Ntrials             = numel(trial_values);
N_CVranges          = size(CV_bounds,1);

% Voltage input to voltage-to-current converter [V]
Vinput_values       = amp_values_actual ./ V2mA;    

subplot_rows        = 3;
subplot_cols        = 1;

% Subplot dimensions to plot all trials for given block in one figure
subplot_rows_all = 3;
subplot_cols_all = 7;

check_Nsubplots = 1;
while (check_Nsubplots == 1)
    if ((subplot_rows_all * subplot_cols_all) < Ntrials)
        subplot_rows_all = subplot_rows_all + 1;
    else
        check_Nsubplots = 0;
    end
end

%% Time and CV parameters
t_bounds_slowest    = stim_rec_dist ./ CV_bounds_slowest;

% Latency bounds
CV_latencybounds    = stim_rec_dist ./ CV_bounds;

stim_period         = 1 / stim_freq;          % [ms]
t_max               = t_bounds_slowest;       % [ms]
earliest_tstart_noise = stim_period - t_max;

pretime             = earliest_tstart_noise;  % Time before event to include in window [ms]
                                                % Negative to specify window starting after event
posttime            = t_max;                  % Time after event to include in window [ms]

% t_window
pretime_idx         = round(pretime * fs / 1000);
posttime_idx        = round(posttime * fs / 1000);
Nidx                = numel(-pretime_idx : posttime_idx);
t_window            = -dt * pretime_idx : dt : dt * posttime_idx;
t_window            = t_window * 1e3;
CV_window           = stim_rec_dist ./ t_window;

% Convert latency bounds to indices
CV_latencybounds_idx = zeros(size(CV_latencybounds));
for CVrange_ind = 1:N_CVranges
    for bound_ind = 1:2
        [~,CV_latencybounds_idx(CVrange_ind,bound_ind)] = min(abs(t_window - CV_latencybounds(CVrange_ind,bound_ind)));
    end
end

% Latency duration for each range of CVs
CV_latencydurations = -diff(CV_bounds, 1, 2);

% Noise for each trial
clear noise_latencybounds
if subtractnoise == 1
    noise_latencybounds(1) = -earliest_tstart_noise;
    noise_latencybounds(2) = -1;
end

noise_durn          = noise_latencybounds(2) - noise_latencybounds(1);
noise_multipliers   = CV_latencydurations ./ noise_durn;

%% Plotting parameters
thresh_targets  = [0.2 0.5 0.8];
cmp_thresh      = parula;
Ncmp            = numel(thresh_targets);
cmp_thresh      = cmp_thresh(1:floor(size(cmp_thresh,1)/Ncmp):end,:);

% Subplot dimensions for a panel for each CV range
subplot_rows_CVranges = 2;
subplot_cols_CVranges = 3;
if ((subplot_rows_CVranges * subplot_cols_CVranges) < N_CVranges)
    error('Need more subplots for N_CVranges')
end

cmp_CVranges    = parula;
Ncmp            = N_CVranges;
cmp_CVranges    = cmp_CVranges(1:floor(size(cmp_CVranges,1)/Ncmp):end,:);

y_CV_bounds = y_min : 1 : y_max;

% Axis limits
CAPvsTime_axlim = [-0.5 max(CV_latencybounds(:)) y_min y_max];
CAPvsCV_axlim   = [0    max(CV_bounds(:))        y_min y_max];
CAPvsFreq_axlim = [0    4000                     0     y_max_f];

%% Load cached raw data (Mac-compatible)
% Load directly from mat_cache_of_raw_data
% Extract the prefix from file_prefix_outputs (e.g., "20220531_1")
[~, prefix_name] = fileparts(file_prefix_outputs);

% Construct path to raw cache
raw_cache_root = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data';
cache_filepath = fullfile(raw_cache_root, prefix_name, sprintf('ch%d.mat', ch_num_neural));

if ~exist(cache_filepath, 'file')
    warning('Cached raw data not found: %s\nRun AnalysisScript_CAPs on Windows first to create the cache.', cache_filepath);
    return;
end

fprintf('Loading cached raw data: %s\n', cache_filepath);
load(cache_filepath, 'ch_Vinput', 'ch_neural', 'cached_trial_values', 'fs', 'gain_value');

% Ensure trial_values match
assert(isequal(trial_values, cached_trial_values), ...
    'Trial values mismatch between TestPulseData and cached data');

% ch_Vinput and ch_neural are already loaded as cell arrays

%% Apply 60Hz filtering to raw trials (filter each trial individually to avoid memory issues)
fprintf('Applying 60Hz filtering to each trial...\n');

% Filter each trial separately
for trial_ind = 1:Ntrials
    neural_raw = ch_neural{trial_ind} * labchart_unit_conv * 1e6; % Convert to uV
    ch_neural{trial_ind} = filter_periodic_artifact(neural_raw, fs, 60);
end

%% Slice epochs from filtered data
fprintf('Slicing epochs from filtered signal...\n');

% Find triggers for each trial
min_samples_between_pulses = round(0.9 * (stim_period/1e3) * fs);
edge_detection_threshold = 0.01;

all_trig_indices = zeros(Ntrig_expected, Ntrials);

for trial_ind = 1:Ntrials
    vinput_trial = ch_Vinput{trial_ind};
    abs_derivative = diff(abs(vinput_trial));
    trig_indices(1) = find(abs_derivative > edge_detection_threshold, 1);
    
    for pulse_idx = 2:Ntrig_expected
        search_start_idx = min_samples_between_pulses + trig_indices(pulse_idx - 1);
        relative_onset = find(abs_derivative(search_start_idx:end) > edge_detection_threshold, 1);
        trig_indices(pulse_idx) = search_start_idx - 1 + relative_onset;
    end
    
    all_trig_indices(:, trial_ind) = trig_indices;
end

% Calculate epoch window
post_samp = posttime_idx;
pre_samp = pretime_idx;
n_epoch = pre_samp + post_samp + 1;

% Find low and high amplitude trials
low_amp_idx = find(amp_values_actual == min(amp_values_actual), 1);
high_amp_idx = find(amp_values_actual == max(amp_values_actual), 1);

fprintf('  Low amp trial: %d (amp=%.3g mA), High amp trial: %d (amp=%.3g mA)\n', ...
    low_amp_idx, amp_values_actual(low_amp_idx), ...
    high_amp_idx, amp_values_actual(high_amp_idx));

%% Build artifact template from low amplitude trial
fprintf('Building artifact template from low-amplitude epochs...\n');

% Slice epochs from low amplitude trial only
trig_offsets = all_trig_indices(:, low_amp_idx);

epochs_low = zeros(n_epoch, Ntrig_expected);
for ep = 1:Ntrig_expected
    s = trig_offsets(ep) - pre_samp;
    e = trig_offsets(ep) + post_samp;
    if s >= 1 && e <= numel(ch_neural{low_amp_idx})
        ep_data = ch_neural{low_amp_idx}(s:e);
        ep_data = ep_data - mean(ep_data); % baseline correct
        epochs_low(:, ep) = ep_data;
    end
end

% Average low-amplitude epochs to create template
low_avg = mean(epochs_low, 2);  % n_epoch x 1

%% Subtract artifact from all trials
fprintf('Subtracting artifact template from all epochs...\n');

% Scale factor for template
amp_scale = amp_values_actual / amp_values_actual(low_amp_idx);

ch_neural_avg = zeros(Nidx, Ntrials);

for trial_ind = 1:Ntrials
    trig_offsets = all_trig_indices(:, trial_ind);
    
    % Slice epochs for this trial
    epochs_trial = zeros(n_epoch, Ntrig_expected);
    for ep = 1:Ntrig_expected
        s = trig_offsets(ep) - pre_samp;
        e = trig_offsets(ep) + post_samp;
        if s >= 1 && e <= numel(ch_neural{trial_ind})
            ep_data = ch_neural{trial_ind}(s:e);
            ep_data = ep_data - mean(ep_data); % baseline correct
            epochs_trial(:, ep) = ep_data;
        end
    end
    
    % Subtract scaled template
    epochs_clean = zeros(n_epoch, Ntrig_expected);
    for ep = 1:Ntrig_expected
        epochs_clean(:, ep) = epochs_trial(:, ep) - amp_scale(trial_ind) * low_avg;
    end
    
    % Average clean epochs
    if CAP_avg_type == 1  % Mean
        ch_neural_avg(:, trial_ind) = mean(epochs_clean, 2);
    elseif CAP_avg_type == 2  % Median
        ch_neural_avg(:, trial_ind) = median(epochs_clean, 2);
    end
end

%% Remove baseline offset
for trial_ind = 1:Ntrials
    y_tmp = ch_neural_avg(:, trial_ind);
    if zerobaseline == 1
        % Baseline window (pre-stimulus)
        baseline_indices = find(t_window < 0);
        if ~isempty(baseline_indices)
            baseline_tmp = mean(y_tmp(baseline_indices));
            ch_neural_avg(:, trial_ind) = y_tmp - baseline_tmp;
        end
    end
end

%% Save output data
% Convert gain from unitless to proper scaling
ch_neural_avg = ch_neural_avg / gain_value;

% Save to artifact-removed output folder
fname = [file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'];
save(fname, 'expdate', 'stim_rec_dist', 'ch_num_neural',...
    'trial_values', 'amp_values_nominal', 'amp_values_actual',...
    't_window', 'CV_window', 'ch_neural_avg', 'subplot_rows_all', 'subplot_cols_all',...
    'CAPvsTime_axlim', 'CAPvsCV_axlim',...
    'CV_bounds', 'CV_latencybounds', 'y_CV_bounds', 't_bounds_slowest',...
    '-v7.3');

fprintf('Saved artifact-removed data to: %s\n', fname);

%% Plot stim-triggered averages vs. t
figure
set(gcf, 'position', scrnsz)

plot_CAPvsTime(Ntrials, ch_neural_avg, subplot_rows_all, subplot_cols_all,...
    t_window, N_CVranges, CV_latencybounds, y_CV_bounds,...
    t_bounds_slowest, CAPvsTime_axlim, expdate, ch_num_neural,...
    phase_type, stim_PW, trial_values, amp_values_nominal)

fname = [file_prefix_outputs '/CAPvsTime_ch' num2str(ch_num_neural)];
saveas(gcf, [fname '.jpg']);
saveas(gcf, [fname '.fig']);

%% Plot stim-triggered averages vs. CV
figure
set(gcf, 'position', scrnsz)

plot_CAPvsCV(Ntrials, ch_neural_avg, subplot_rows_all, subplot_cols_all,...
    CV_window, N_CVranges, CV_bounds, y_CV_bounds,...
    CV_bounds_slowest, CAPvsCV_axlim, expdate, ch_num_neural,...
    phase_type, stim_PW, trial_values, amp_values_nominal)

fname = [file_prefix_outputs '/CAPvsCV_ch' num2str(ch_num_neural)];
saveas(gcf, [fname '.jpg']);
saveas(gcf, [fname '.fig']);

end
