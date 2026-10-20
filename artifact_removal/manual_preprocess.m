clear
clc
close all
% manual_preprocess.m
% Load a preprocessed CAP .mat file and the corresponding raw data cache
% .mat file, inspect fields, and plot raw vs preprocessed signals.

% preproc_path = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data/20221201_2/CAPData_ch1.mat';
% raw_path     = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20221201_2/ch1.mat';
preproc_path = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data/20220518_1/CAPData_ch1.mat';
raw_path     = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220518_1/ch1.mat';

% ── Load ─────────────────────────────────────────────────────────────────────
preproc = load(preproc_path);
raw     = load(raw_path);

% ── Print fields ─────────────────────────────────────────────────────────────
print_fields(preproc, preproc_path);
print_fields(raw,     raw_path);

% ── Plot: single amplitude level ─────────────────────────────────────────────
% Change AMP_IDX to inspect a different amplitude level (1 = lowest, n_amp = highest)
[~, AMP_IDX] = min(preproc.amp_values_actual); %plot the max amplitude
%AMP_IDX = 3;
n_amp = numel(raw.ch_neural);
t_win = preproc.t_window;
fs    = raw.fs;

% Find stimulus onset via Vinput threshold
vinput_ref  = raw.ch_Vinput{1};
thresh      = 0.01 * max(abs(vinput_ref(:)));
stim_sample = find(abs(vinput_ref(:, 1)) > thresh, 1, 'first');
stim_ms     = (stim_sample - 1) / fs * 1000;
fprintf('Stimulus onset detected at %.2f ms into raw sweep\n', stim_ms);

% Crop raw to same window as t_window
t_start_ms = t_win(1);
t_end_ms   = t_win(end);
crop_start  = stim_sample + round(t_start_ms / 1000 * fs);
crop_end    = stim_sample + round(t_end_ms   / 1000 * fs);
t_raw_crop  = (0 : crop_end - crop_start) / fs * 1000 + t_start_ms;

ZOOM_XLIM  = [-0.5, 2.0];   % ms
amp_label  = sprintf('%.3g mA', preproc.amp_values_actual(AMP_IDX));
raw_full   = raw.ch_neural{AMP_IDX};    % full sweep for this amplitude
vinput_amp = raw.ch_Vinput{AMP_IDX};    % Vinput for this amplitude
Ntrig_expected = 10;

% Find all 10 trigger times in this amplitude's Vinput
abs_deriv  = diff(abs(vinput_amp));
edge_thresh = 0.01;
min_gap     = round(0.9 * (1/0.003096) * fs / 1000);  % 90% of stim period in samples

trig_indices = zeros(1, Ntrig_expected);
trig_indices(1) = find(abs_deriv > edge_thresh, 1, 'first');
for p = 2:Ntrig_expected
    search_start = trig_indices(p-1) + min_gap;
    rel = find(abs_deriv(search_start:end) > edge_thresh, 1, 'first');
    trig_indices(p) = search_start - 1 + rel;
end
fprintf('Found %d triggers for amp %s\n', Ntrig_expected, amp_label);

% Slice each epoch around its own trigger
n_samp_epoch = crop_end - crop_start + 1;
epochs = zeros(n_samp_epoch, Ntrig_expected);
for p = 1:Ntrig_expected
    s = trig_indices(p) + round(t_start_ms / 1000 * fs);
    e = s + n_samp_epoch - 1;
    if s < 1 || e > numel(raw_full)
        epochs(:, p) = NaN;
    else
        epochs(:, p) = raw_full(s:e);
    end
end

% Convert to µV (same as pipeline: multiply by labchart_unit_conv=1, then ×1e6, divide by gain)
epochs_uV = epochs * 1e6 / raw.gain_value;

% ── Figure 1: 10 sliced epochs (individual trials) ───────────────────────────
figure('Name', sprintf('10 Epochs — %s', amp_label), ...
       'NumberTitle', 'off', 'Position', [50 100 1200 600]);

cmap_trials = parula(Ntrig_expected);
ax_epochs = gobjects(1, Ntrig_expected);
for p = 1:Ntrig_expected
    ax_epochs(p) = subplot(2, 5, p);
    plot(t_raw_crop, epochs_uV(:, p), 'Color', cmap_trials(p,:), 'LineWidth', 1.2);
    title(sprintf('Trial %d', p), 'FontSize', 10);
    xlim(ZOOM_XLIM);
    grid on;
    if p == 1 || p == 6
        ylabel('µV', 'FontSize', 10);
    end
    if p > 5
        xlabel('Time (ms)', 'FontSize', 10);
    end
end
linkaxes(ax_epochs, 'xy');
sgtitle(sprintf('Individual epochs — %s, ch%d, dist=%g mm', ...
    amp_label, preproc.ch_num_neural, preproc.stim_rec_dist), 'FontSize', 13);

% ── Figure 2: all epochs overlaid + preprocessed median ──────────────────────
figure('Name', sprintf('Raw vs Preprocessed — %s', amp_label), ...
       'NumberTitle', 'off', 'Position', [100 100 900 600]);

ax1 = subplot(2, 1, 1);
hold on;
for p = 1:Ntrig_expected
    plot(t_raw_crop, epochs_uV(:, p), 'Color', [0.2 0.5 0.9 0.4], 'LineWidth', 0.8);
end
hold off;
title(sprintf('All 10 epochs overlaid — %s', amp_label), 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 11);
xlim(ZOOM_XLIM);
grid on;

ax2 = subplot(2, 1, 2);
plot(t_win, preproc.ch_neural_avg(:, AMP_IDX), 'Color', [0.85 0.2 0.2], 'LineWidth', 1.8);
title(sprintf('Preprocessed median — %s', amp_label), 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 11);
xlabel('Time (ms)', 'FontSize', 11);
xlim(ZOOM_XLIM);
grid on;

linkaxes([ax1, ax2], 'x');
sgtitle(sprintf('ch%d, dist=%g mm', preproc.ch_num_neural, preproc.stim_rec_dist), ...
        'FontSize', 13);

% for amp_idx = 1:n_amp
%     amp_label  = sprintf('%.3g mA', preproc.amp_values_actual(amp_idx));
%     col        = cmap(amp_idx, :);
%     raw_trials = raw.ch_neural{amp_idx};
%     n_trials   = size(raw_trials, 2);
% 
%     ax_raw(amp_idx) = subplot(2, n_amp, amp_idx);
%     hold on;
%     for tr = 1:n_trials
%         seg = raw_trials(crop_start:crop_end, tr);
%         plot(t_raw_crop, seg, 'Color', [col, 0.4], 'LineWidth', 0.8);
%     end
%     hold off;
%     title(amp_label, 'FontSize', 10);
%     grid on;
%     if amp_idx == 1; ylabel('Raw (µV)', 'FontSize', 11); end
%     %xlim(ZOOM_XLIM);
% 
%     ax_preproc(amp_idx) = subplot(2, n_amp, n_amp + amp_idx);
%     plot(t_win, preproc.ch_neural_avg(:, amp_idx), 'Color', col, 'LineWidth', 1.5);
%     grid on;
%     xlabel('Time (ms)', 'FontSize', 10);
%     if amp_idx == 1; ylabel('Preprocessed (µV)', 'FontSize', 11); end
%     %xlim(ZOOM_XLIM);
% end
% linkaxes([ax_raw, ax_preproc], 'x');


% ── Helper ───────────────────────────────────────────────────────────────────
function print_fields(data, label)
    fprintf('\n=== Fields in %s ===\n', label);
    fields = fieldnames(data);
    for i = 1:numel(fields)
        f   = fields{i};
        val = data.(f);
        sz  = size(val);
        sz_str = strjoin(arrayfun(@num2str, sz, 'UniformOutput', false), ' x ');

        if isnumeric(val) && numel(val) == 1
            fprintf('  %-30s [%s]  =  %g\n', f, sz_str, val);
        elseif isnumeric(val) && numel(val) <= 10
            fprintf('  %-30s [%s]  =  %s\n', f, sz_str, mat2str(val, 4));
        elseif ischar(val) || (isstring(val) && numel(val) == 1)
            fprintf('  %-30s [%s]  =  "%s"\n', f, sz_str, val);
        else
            fprintf('  %-30s [%s]\n', f, sz_str);
        end
    end
end