% visualize_amplitudes.m
%
% Loads raw cache for a given dataset folder + channel, slices the 10
% triggered epochs at every stimulation amplitude, averages them, and
% plots all averaged CAPs overlaid and as a stacked grid so you can
% inspect how the waveform changes with increasing stimulus amplitude.
%
% ── User settings ─────────────────────────────────────────────────────────────
DATASET_FOLDER = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20221122_1';
CH_NUM         = 1;
NTRIG_EXPECTED = 10;
ZOOM_XLIM      = [-0.5, 3.0];   % ms
% ──────────────────────────────────────────────────────────────────────────────

%% Load
raw_path = fullfile(DATASET_FOLDER, sprintf('ch%d.mat', CH_NUM));
fprintf('Loading: %s\n', raw_path);
raw = load(raw_path);

fs           = raw.fs;
dt_ms        = 1000 / fs;
amp_nominal  = raw.amp_values_nominal(:)';
n_amp        = numel(raw.ch_neural);

fprintf('n_amp=%d  amps: %s mA\n', n_amp, num2str(amp_nominal, '%.3g '));

%% Detect triggers and slice epochs
stim_period_ms = 1000 / 3.096;
min_gap        = round(0.9 * stim_period_ms / 1000 * fs);
edge_thresh    = 0.01;

pretime_ms  = 1.0;
posttime_ms = 5.0;
pre_samp    = round(pretime_ms  / dt_ms);
post_samp   = round(posttime_ms / dt_ms);
n_epoch     = pre_samp + post_samp + 1;
t_epoch     = (-pre_samp:post_samp) * dt_ms;

avg_epochs = zeros(n_epoch, n_amp);   % one averaged CAP per amplitude

for ai = 1:n_amp
    neural = raw.ch_neural{ai};
    vinput = raw.ch_Vinput{ai};
    trigs  = find_triggers(vinput, NTRIG_EXPECTED, min_gap, edge_thresh);

    ep = zeros(n_epoch, NTRIG_EXPECTED);
    for ti = 1:NTRIG_EXPECTED
        s = trigs(ti) - pre_samp;
        e = trigs(ti) + post_samp;
        if s >= 1 && e <= numel(neural)
            trial = neural(s:e) * 1e6 / raw.gain_value;
            trial = trial - mean(trial(t_epoch < 0));   % baseline correct
            ep(:, ti) = trial;
        end
    end
    avg_epochs(:, ai) = mean(ep, 2);
end

%% Plot 1: all amplitudes overlaid, colored by amplitude
zoom_mask = t_epoch >= ZOOM_XLIM(1) & t_epoch <= ZOOM_XLIM(2);
cmap      = parula(n_amp);

figure('Name', sprintf('All amplitudes overlaid — %s ch%d', ...
    DATASET_FOLDER(end-10:end), CH_NUM), ...
    'NumberTitle', 'off', 'Position', [50 80 850 500]);
hold on;
for ai = 1:n_amp
    plot(t_epoch(zoom_mask), avg_epochs(zoom_mask, ai), ...
         'Color', cmap(ai,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('%.3g mA', amp_nominal(ai)));
end
xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'HandleVisibility', 'off');
hold off;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Amplitude (µV)', 'FontSize', 12);
title(sprintf('Averaged CAPs — all amplitudes  ch%d', CH_NUM), 'FontSize', 13);
legend('Location', 'eastoutside', 'FontSize', 8);

sm = colorbar; sm.Label.String = 'Stim amplitude (mA)';
colormap(parula); clim([amp_nominal(1) amp_nominal(end)]);
grid on;

%% Plot 2: stacked grid — one subplot per amplitude
ncols = min(n_amp, 5);
nrows = ceil(n_amp / ncols);

% Global y limits across all amplitudes
y_all = avg_epochs(zoom_mask, :);
ylims = [min(y_all(:)) - 0.05*range(y_all(:)), ...
         max(y_all(:)) + 0.05*range(y_all(:))];

figure('Name', sprintf('Stacked grid — %s ch%d', ...
    DATASET_FOLDER(end-10:end), CH_NUM), ...
    'NumberTitle', 'off', 'Position', [100 80 250*ncols 200*nrows]);

ax_all = gobjects(nrows, ncols);
for ai = 1:n_amp
    r = ceil(ai / ncols);
    c = mod(ai - 1, ncols) + 1;
    ax = subplot(nrows, ncols, ai);
    ax_all(r, c) = ax;

    plot(t_epoch(zoom_mask), avg_epochs(zoom_mask, ai), ...
         'Color', cmap(ai,:), 'LineWidth', 1.3);
    hold on;
    xline(0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
    hold off;
    title(sprintf('%.3g mA', amp_nominal(ai)), 'FontSize', 9);
    xlim(ZOOM_XLIM);
    ylim(ylims);
    grid on;
    if c == 1; ylabel('µV', 'FontSize', 8); end
    if r == nrows; xlabel('Time (ms)', 'FontSize', 8); end
end

% Remove empty subplots
for idx = n_amp+1 : nrows*ncols
    r = ceil(idx / ncols);
    c = mod(idx - 1, ncols) + 1;
    if ishandle(ax_all(r,c)); delete(ax_all(r,c)); end
end

sgtitle(sprintf('Averaged CAPs per amplitude — ch%d', CH_NUM), 'FontSize', 13);

%% Local functions
function trig = find_triggers(vinput_col, ntrig, min_gap, edge_thr)
    abs_d   = diff(abs(vinput_col));
    trig    = zeros(1, ntrig);
    trig(1) = find(abs_d > edge_thr, 1, 'first');
    for p   = 2:ntrig
        s   = trig(p-1) + min_gap;
        r   = find(abs_d(s:end) > edge_thr, 1, 'first');
        trig(p) = s - 1 + r;
    end
end
