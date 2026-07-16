% artifact_removal.m
%
% Loads all preprocessed CAPs for a given expdate and nerve type,
% aligns the stimulus artifacts across channels, normalises amplitudes,
% fits a template artifact (biphasic/monophasic square wave + exponential
% decay tail) to the artifact region, subtracts it, and plots before/after.
%
% stim_PW and phase_type are read automatically from the xlsx metadata.
%
% --- USER SETTINGS -----------------------------------------------------------
TARGET_EXPDATE  = '20221201';
TARGET_NERVE    = 'SN';          % 'SN' | 'cVN' | 'tibial'
AMP_IDX         = [];            % [] = max amplitude; or set e.g. 5

PREPROC_ROOT    = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data';
ZOOM_XLIM       = [-0.5, 3.0];   % ms — plot window around stimulus
% -----------------------------------------------------------------------------

%% 0. Read stim_PW and phase_type from xlsx
% Parse raw_data_path from local_config.yml
local_config_file = 'local_config.yml';
fid = fopen(local_config_file);
raw_data_path = '';
line_i = fgetl(fid);
while ischar(line_i)
    if contains(line_i, 'raw_data_path: ')
        raw_data_path = strtrim(line_i(length('raw_data_path: ')+2:end-1));
    end
    line_i = fgetl(fid);
end
fclose(fid);

XLSX_PATH = fullfile(raw_data_path, 'AnalysisScript_CAPs_Batch_parameters.xlsx');
meta = readtable(XLSX_PATH);
meta.expdate_str = arrayfun(@(x) num2str(x), meta.expdate_values, ...
                             'UniformOutput', false);

meta_match = meta(strcmp(meta.expdate_str, TARGET_EXPDATE) & ...
                  strcmp(meta.rec_nerve,    TARGET_NERVE),   :);

if isempty(meta_match)
    error('No metadata rows found for expdate=%s nerve=%s', ...
          TARGET_EXPDATE, TARGET_NERVE);
end

phase_types = unique(meta_match.phase_type);
pw_values   = unique(meta_match.stim_PW);

fprintf('\nMetadata for %s %s:\n', TARGET_EXPDATE, TARGET_NERVE);
fprintf('  phase_type(s): %s\n', strjoin(phase_types, ', '));
fprintf('  stim_PW(s):    %s ms\n', num2str(pw_values(:)', '%.3f '));

% Filter out biphasic rows — keep monophasic only for artifact template
mono_mask   = strcmp(meta_match.phase_type, 'monophasic');
if any(mono_mask)
    meta_match = meta_match(mono_mask, :);
    fprintf('Filtered to monophasic only (%d rows remaining)\n', sum(mono_mask));
end

% Recompute after filtering
phase_types = unique(meta_match.phase_type);
pw_values   = unique(meta_match.stim_PW);

PHASE_TYPE = phase_types{1};
STIM_PW_MS = pw_values(1);
fprintf('Using phase_type="%s", stim_PW=%.3f ms\n\n', PHASE_TYPE, STIM_PW_MS);

%% 1. Discover CAPData_ch*.mat files for this expdate + nerve only
% The xlsx already tells us which folder prefixes belong to TARGET_NERVE —
% use those to whitelist folders before loading anything.
valid_prefixes = meta_match.file_prefix_outputs_values;   % e.g. {'20221201_7','20221201_8',...}

all_files = [];
for pi = 1:numel(valid_prefixes)
    prefix  = valid_prefixes{pi};
    pattern = fullfile(PREPROC_ROOT, prefix, 'CAPData_ch*.mat');
    found   = dir(pattern);
    all_files = [all_files; found];  %#ok<AGROW>
end

if isempty(all_files)
    error('No CAPData files found for expdate=%s nerve=%s', ...
          TARGET_EXPDATE, TARGET_NERVE);
end
fprintf('Found %d mat file(s) for %s %s\n', numel(all_files), TARGET_EXPDATE, TARGET_NERVE);

%% 2. Load data from each file
records = struct([]);
for fi = 1:numel(all_files)
    fpath = fullfile(all_files(fi).folder, all_files(fi).name);
    d = load(fpath, 'ch_neural_avg', 't_window', 'amp_values_actual', ...
                    'ch_num_neural', 'stim_rec_dist', 'expdate');

    % Pick amplitude column
    if isempty(AMP_IDX)
        [~, ai] = max(d.amp_values_actual);
    else
        ai = AMP_IDX;
    end

    % Look up dataset_stim_polarity for this folder from the xlsx
    % all_files(fi).folder is the full path ending in the prefix, e.g. '.../20221201_3'
    [~, prefix] = fileparts(all_files(fi).folder);   % extracts '20221201_3'
    row_idx = find(strcmp(meta_match.file_prefix_outputs_values, prefix), 1);
    if ~isempty(row_idx)
        polarity = meta_match.dataset_stim_polarity(row_idx);
        if iscell(polarity); polarity = polarity{1}; end
        polarity = double(polarity);
        if isnan(polarity); polarity = 1; end
    else
        polarity = 1;
    end

    signal = d.ch_neural_avg(:, ai);

    % Negate polarity 2 signals so all recordings face the same direction
    if polarity == 2
        signal = -signal;
        fprintf('  Flipping polarity 2: %s / %s\n', prefix, all_files(fi).name);
    end

    records(fi).signal       = signal;
    records(fi).t_window     = d.t_window(:);
    records(fi).ch           = d.ch_num_neural;
    records(fi).dist_mm      = d.stim_rec_dist;
    records(fi).amp_mA       = d.amp_values_actual(ai);
    records(fi).polarity      = polarity;
    records(fi).fname        = all_files(fi).name;
    records(fi).folder       = all_files(fi).folder;
end

N = numel(records);
fprintf('Loaded %d channels (%d polarity-2 flipped)\n', N, ...
    sum([records.polarity] == 2));

%% 3. Artifact region definition
% Artifact spans from t=0 to end of pulse(s) + ~0.5 ms ringdown buffer
if strcmp(PHASE_TYPE, 'biphasic')
    artifact_end_ms = STIM_PW_MS * 2 + 0.5;
else
    artifact_end_ms = STIM_PW_MS + 0.5;
end

t_ref      = records(1).t_window;
art_mask   = (t_ref >= 0) & (t_ref <= artifact_end_ms);
pre_mask   = (t_ref >= -1.0) & (t_ref < 0);   % 1 ms pre-stim baseline

fprintf('artifact_end_ms=%.3f  art_mask samples=%d  pre_mask samples=%d\n', ...
    artifact_end_ms, sum(art_mask), sum(pre_mask));

% Safety check
if sum(art_mask) < 5
    error('art_mask has only %d samples — check STIM_PW_MS (%.3f) and PHASE_TYPE (%s)', ...
        sum(art_mask), STIM_PW_MS, PHASE_TYPE);
end

% Quick sanity: print peak value of first record's raw signal
fprintf('Record 1 raw signal range: [%.2f, %.2f] uV\n', ...
    min(records(1).signal), max(records(1).signal));

%% 4. Baseline-correct and amplitude-normalise each channel
signals_norm = zeros(numel(t_ref), N);
peak_arts    = zeros(1, N);

for i = 1:N
    sig = records(i).signal;
    % Subtract pre-stimulus mean
    sig = sig - mean(sig(pre_mask));
    % Normalise by peak absolute value in artifact window
    peak_arts(i) = max(abs(sig(art_mask)));
    if peak_arts(i) > 0
        signals_norm(:, i) = sig ./ peak_arts(i);
    else
        signals_norm(:, i) = sig;
    end
end

%% 5. Build artifact template
% Strategy: the monophasic artifact is a square pulse (0 to stim_PW) 
% followed by an RC-style exponential decay. We know the pulse width, so:
%   - During the pulse:  use the consensus average of the normalised signals
%   - After the pulse:   fit A * exp(-t/tau) with t=0 at stim_PW
% This requires fitting only tau (one parameter).

template_avg = mean(signals_norm(art_mask, :), 2);
t_art        = t_ref(art_mask)';

% Pulse end index within the artifact window
pw_idx = find(t_art >= STIM_PW_MS, 1, 'first');
if isempty(pw_idx); pw_idx = numel(t_art); end

% Post-pulse tail: everything from pulse end to artifact window end
t_tail   = t_art(pw_idx:end) - t_art(pw_idx);   % zero-based from pulse end
y_tail   = template_avg(pw_idx:end);
A0       = mean(y_tail(1:min(3,end)));           % amplitude at pulse end

% Fit tau via log-linear regression on abs(y_tail)
sign_tail   = sign(A0);
y_abs       = abs(y_tail);
noise_floor = 0.01 * max(abs(template_avg));
valid       = y_abs > noise_floor;

if sum(valid) >= 2
    p       = polyfit(t_tail(valid), log(y_abs(valid)), 1);
    tau_fit = -1 / p(1);
    if tau_fit <= 0 || tau_fit > 50
        tau_fit = STIM_PW_MS * 3;   % fallback: decay over ~3x the pulse width
    end
else
    tau_fit = STIM_PW_MS * 3;
end

fprintf('Pulse width=%.3f ms  Decay tau=%.3f ms\n', STIM_PW_MS, tau_fit);

% Assemble full template over the artifact window:
%   pulse body = data average, post-pulse = exponential from A0
template_full = zeros(size(t_ref));
template_full(art_mask) = template_avg;   % use real data for pulse body

% Overwrite the post-pulse portion with the fitted exponential
post_pulse_in_art = art_mask & (t_ref >= STIM_PW_MS);
t_pp = t_ref(post_pulse_in_art)' - STIM_PW_MS;
template_full(post_pulse_in_art) = A0 * exp(-t_pp ./ tau_fit);

% Extend decay beyond the artifact window
post_mask = t_ref > artifact_end_ms;
t_post    = t_ref(post_mask)' - STIM_PW_MS;
template_full(post_mask) = A0 * exp(-t_post ./ tau_fit);

%% 6. Subtract scaled template, then restore original scale and polarity
signals_clean = zeros(size(signals_norm));
for i = 1:N
    % Work in normalised space: subtract normalised template
    signals_clean(:, i) = signals_norm(:, i) - template_full;

    % Rescale back to original µV amplitude
    signals_clean(:, i) = signals_clean(:, i) * peak_arts(i);

    % Restore original polarity (undo the polarity-2 flip)
    if records(i).polarity == 2
        signals_clean(:, i) = -signals_clean(:, i);
    end
end

%% 7. Plot: one figure per channel — before (with template) on top, after on bottom
cmap     = lines(N);
zoom_idx = (t_ref >= ZOOM_XLIM(1)) & (t_ref <= ZOOM_XLIM(2));

for i = 1:N
    fig_title = sprintf('ch%d  %.1f mm  %.3g mA  pol%d — %s %s', ...
        records(i).ch, records(i).dist_mm, records(i).amp_mA, ...
        records(i).polarity, TARGET_EXPDATE, TARGET_NERVE);

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [50 + (i-1)*30, 80, 700, 550]);

    % Top: before + template overlay (show original polarity)
    ax1 = subplot(2, 1, 1);
    sig_bc = records(i).signal - mean(records(i).signal(pre_mask));
    % Restore original polarity for display
    if records(i).polarity == 2
        sig_bc = -sig_bc;
    end
    plot(t_ref(zoom_idx), sig_bc(zoom_idx), ...
         'Color', cmap(i,:), 'LineWidth', 1.4);
    hold on;
    % Template scaled back and shown in original polarity
    template_scaled = template_full * peak_arts(i);
    if records(i).polarity == 2; template_scaled = -template_scaled; end
    plot(t_ref(zoom_idx), template_scaled(zoom_idx), ...
         'k--', 'LineWidth', 1.2);
    xline(0,               'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    xline(artifact_end_ms, 'Color', [0.8 0.3 0.3], 'LineStyle', ':');
    hold off;
    ylabel('Before (µV)',   'FontSize', 11);
    legend({'Signal', 'Template'}, 'FontSize', 9, 'Location', 'northeast');
    title(fig_title, 'FontSize', 11);
    grid on;
    xlim(ZOOM_XLIM);

    % Bottom: after subtraction
    ax2 = subplot(2, 1, 2);
    plot(t_ref(zoom_idx), signals_clean(zoom_idx, i), ...
         'Color', cmap(i,:), 'LineWidth', 1.4);
    hold on;
    xline(0,               'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    xline(artifact_end_ms, 'Color', [0.8 0.3 0.3], 'LineStyle', ':');
    hold off;
    ylabel('After subtraction (µV)', 'FontSize', 11);
    xlabel('Time (ms)',              'FontSize', 11);
    grid on;
    xlim(ZOOM_XLIM);

    linkaxes([ax1, ax2], 'x');
end

%% 8. Plot template alone
figure('Name', 'Artifact Template', 'NumberTitle', 'off', 'Position', [100 700 700 300]);
plot(t_ref(zoom_idx), template_full(zoom_idx), 'k-', 'LineWidth', 1.8);
xline(0,              'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'Label', 't=0');
xline(artifact_end_ms,'Color', [0.8 0.3 0.3], 'LineStyle', ':', 'Label', 'artifact end');
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Normalised amplitude', 'FontSize', 12);
title(sprintf('Consensus artifact template  (tau = %.2f ms)', tau_fit), 'FontSize', 13);
grid on;
xlim(ZOOM_XLIM);
