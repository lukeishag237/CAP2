% fit_artifact.m
%
% Artifact subtraction using raw data only (no preprocessed mat needed).
%
% Strategy:
%   1. Load ch1 raw data from the mat_cache for a given dataset folder
%   2. Detect stimulus trigger times from ch_Vinput
%   3. At the LOWEST amplitude (assume pure artifact, no CAP): slice the
%      10 triggered epochs and average them → artifact template
%   4. Scale the template to match the artifact amplitude at the HIGHEST
%      amplitude, then subtract → cleaned signal
%   5. Plot: individual low-amp epochs | avg template | high-amp before & after

% ── User settings ─────────────────────────────────────────────────────────────
addpath('/Users/lukeishag/Documents/PERL/CAP2/external_dependencies/custom_denoising_scripts');
%DATASET_FOLDER = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20221201_2';
DATASET_FOLDER = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220518_6';


CH_NUM         = 2;         % which channel to use (ch1 = ch_neural{} from ch1.mat)
NTRIG_EXPECTED = 10;        % pulses per amplitude level
ZOOM_XLIM      = [-0.5, 3.0]; % ms around stimulus onset for plots
% ──────────────────────────────────────────────────────────────────────────────

%% 1. Load raw cache
raw_path = fullfile(DATASET_FOLDER, sprintf('ch%d.mat', CH_NUM));
fprintf('Loading: %s\n', raw_path);
raw = load(raw_path);

fs              = raw.fs;                     % samples/s (100 000)
dt_ms           = 1000 / fs;                  % ms per sample
amp_nominal     = raw.amp_values_nominal(:)'; % 1 x N_amp
n_amp           = numel(raw.ch_neural);

fprintf('fs=%g Hz  n_amp=%d  amps: %s mA\n', fs, n_amp, num2str(amp_nominal, '%.3g '));

%% 2. Find stimulus onset in Vinput for all amplitudes
% Use the same edge-detection logic as AnalysisScript_CAPs.m:
%   abs_derivative threshold on ch_Vinput, with minimum-distance constraint.

stim_period_ms  = 1 / 0.003096;              % ~323 ms (3.096 Hz stim)
min_gap_samp    = round(0.9 * stim_period_ms / 1000 * fs);
edge_thresh     = 0.01;




%% 3. Slice epochs for every amplitude, channel 1 only
% Epoch window: same approach as manual_preprocess — use the first Vinput
% to determine pretime, then apply that window to all amplitudes.

% Determine pretime from the first trigger of the first amplitude
vinput_ref  = raw.ch_Vinput{1};
trig_ref    = find_triggers(vinput_ref, 1, min_gap_samp, edge_thresh);
stim_samp   = trig_ref(1);
stim_ms     = (stim_samp - 1) / fs * 1000;

% Use a fixed pre/post window (same as manual_preprocess)
% Derive from the preprocessed t_window if available, otherwise use defaults
pretime_ms  = 1.0;   % ms before stimulus
posttime_ms = 5.0;   % ms after stimulus (enough to see CAP for fast fibers)

pre_samp    = round(pretime_ms  / 1000 * fs);
post_samp   = round(posttime_ms / 1000 * fs);
n_epoch     = pre_samp + post_samp + 1;
t_epoch     = (-pre_samp:post_samp) * dt_ms;   % ms, t=0 at stimulus

fprintf('Epoch window: [%.1f, %.1f] ms  (%d samples)\n', ...
    t_epoch(1), t_epoch(end), n_epoch);

% Slice epochs for each amplitude
epochs_all = cell(1, n_amp);   % each cell: n_epoch x NTRIG_EXPECTED

for ai = 1:n_amp
    neural  = raw.ch_neural{ai};
    time = ((1:length(neural))-1)*(1/fs);
    [filtered_neural, time_vector] = filter_periodic_artifact(neural, fs, 60);
    vinput  = raw.ch_Vinput{ai};
    trigs   = find_triggers(vinput, NTRIG_EXPECTED, min_gap_samp, edge_thresh);

    ep = zeros(n_epoch, NTRIG_EXPECTED);
    for ti = 1:NTRIG_EXPECTED
        s = trigs(ti) - pre_samp;
        e = trigs(ti) + post_samp;
        if s < 1 || e > numel(filtered_neural)
            ep(:, ti) = NaN;
        else
            ep(:, ti) = filtered_neural(s:e);
        end
    end
    % Convert to µV: ×1e6 for V→µV, ÷ gain_value to remove amplifier gain
    epochs_all{ai} = ep * 1e6 / raw.gain_value;
end

%% 4. Build artifact template from lowest amplitude (pure artifact)
[~, low_ai]  = min(amp_nominal);
%low_ai  = 2;
[~, high_ai] = max(amp_nominal);

low_epochs   = epochs_all{low_ai};    % n_epoch x 10
template     = mean(low_epochs, 2);   % average across 10 trials → n_epoch x 1

fprintf('Template from amp %.3g mA (index %d)\n', amp_nominal(low_ai), low_ai);
fprintf('Subtract from    amp %.3g mA (index %d)\n', amp_nominal(high_ai), high_ai);

%% 5. Scale template by amplitude ratio, then subtract
high_epochs = epochs_all{high_ai};

% Scale factor = ratio of stimulation amplitudes (stored in mat file).
% Artifact amplitude scales linearly with stim current, so this is the
% physically correct scaling — no peak-fitting needed.
amp_scale = amp_nominal(high_ai) / amp_nominal(low_ai);
fprintf('Amplitude scale factor: %.3g / %.3g = %.4f\n', ...
    amp_nominal(high_ai), amp_nominal(low_ai), amp_scale);

high_clean = zeros(n_epoch, NTRIG_EXPECTED);
for ti = 1:NTRIG_EXPECTED
    high_clean(:, ti) = high_epochs(:, ti) - amp_scale * template;
end

high_clean_avg = mean(high_clean, 2);
high_raw_avg   = mean(high_epochs, 2);

%% 6. Plots

zoom_mask = t_epoch >= ZOOM_XLIM(1) & t_epoch <= ZOOM_XLIM(2);
cmap_gray = @(n) repmat(linspace(0.7, 0.1, n)', 1, 3);   % dark to light gray
col_template = [0.85 0.2 0.2];   % red
col_clean    = [0.1  0.5 0.9];   % blue

% ── Figure 1: Low-amplitude epochs + template ─────────────────────────────────
figure('Name', sprintf('Artifact template — %.3g mA', amp_nominal(low_ai)), ...
       'NumberTitle', 'off', 'Position', [50 100 800 500]);
hold on;
gc = cmap_gray(NTRIG_EXPECTED);
for ti = 1:NTRIG_EXPECTED
    plot(t_epoch(zoom_mask), low_epochs(zoom_mask, ti), ...
         'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
end
plot(t_epoch(zoom_mask), template(zoom_mask), ...
     'Color', col_template, 'LineWidth', 2.0, 'DisplayName', 'Mean template');
xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
hold off;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 12);
title(sprintf('Low-amp epochs + template  (%.3g mA, ch%d)', ...
    amp_nominal(low_ai), CH_NUM), 'FontSize', 13);
legend({'Trials (gray)', 'Mean template'}, 'FontSize', 10);
xlim(ZOOM_XLIM); grid on;

% ── Figure 2: High-amplitude — individual trials before & after ───────────────
figure('Name', sprintf('High-amp trials — %.3g mA', amp_nominal(high_ai)), ...
       'NumberTitle', 'off', 'Position', [100 100 1100 500]);

subplot(1, 2, 1);
hold on;
for ti = 1:NTRIG_EXPECTED
    plot(t_epoch(zoom_mask), high_epochs(zoom_mask, ti), ...
         'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
end
plot(t_epoch(zoom_mask), high_raw_avg(zoom_mask), ...
     'b-', 'LineWidth', 2, 'DisplayName', 'Mean');
xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
hold off;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 12);
title('Before subtraction', 'FontSize', 13);
xlim(ZOOM_XLIM); grid on;

subplot(1, 2, 2);
hold on;
for ti = 1:NTRIG_EXPECTED
    plot(t_epoch(zoom_mask), high_clean(zoom_mask, ti), ...
         'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
end
plot(t_epoch(zoom_mask), high_clean_avg(zoom_mask), ...
     'Color', col_clean, 'LineWidth', 2, 'DisplayName', 'Mean');
xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
hold off;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 12);
title('After template subtraction', 'FontSize', 13);
xlim(ZOOM_XLIM); grid on;

sgtitle(sprintf('%.3g mA  ch%d  — artifact template subtraction', ...
    amp_nominal(high_ai), CH_NUM), 'FontSize', 14);

% ── Figure 3: Mean before vs mean after overlay ───────────────────────────────
figure('Name', 'Before vs After (means)', 'NumberTitle', 'off', ...
       'Position', [150 150 750 450]);
hold on;
plot(t_epoch(zoom_mask), high_raw_avg(zoom_mask),   'k-',  'LineWidth', 1.8, ...
     'DisplayName', sprintf('Before (%.3g mA)', amp_nominal(high_ai)));
plot(t_epoch(zoom_mask), high_clean_avg(zoom_mask), ...
     'Color', col_clean, 'LineWidth', 1.8, ...
     'DisplayName', 'After subtraction');
xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
hold off;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 12);
title(sprintf('Mean CAP before vs after  ch%d', CH_NUM), 'FontSize', 13);
legend('FontSize', 11);
xlim(ZOOM_XLIM); grid on;

%% Local functions
function trig = find_triggers(vinput_col, ntrig, min_gap, edge_thr)
    abs_d = diff(abs(vinput_col));
    trig  = zeros(1, ntrig);
    trig(1) = find(abs_d > edge_thr, 1, 'first');
    for p = 2:ntrig
        s = trig(p-1) + min_gap;
        r = find(abs_d(s:end) > edge_thr, 1, 'first');
        trig(p) = s - 1 + r;
    end
end
