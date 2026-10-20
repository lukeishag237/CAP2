% fit_artifact_batch.m
%
% Batch artifact subtraction over all contacts (ch1..ch3, whichever exist)
% and both stimulation polarities.
%
% For each (polarity, contact):
%   1. Load chN.mat from that polarity's dataset folder
%   2. Detect triggers from ch_Vinput, slice epochs for every amplitude
%   3. Template = mean of the lowest-amplitude epochs (assumed pure artifact)
%   4. Scale template by amplitude ratio, subtract from highest-amplitude epochs
%   5. One figure per (polarity, contact) + one summary figure
%
% Results are collected in `results{pol, ch}` and saved to a .mat file.

% ── User settings ─────────────────────────────────────────────────────────────
addpath('/Users/lukeishag/Documents/PERL/CAP2/external_dependencies/custom_denoising_scripts');

POLARITY_FOLDERS = { ...
    '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220518_5', ... % polarity 1
    '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220518_6'  ... % polarity 2  <-- EDIT
    };
POLARITY_LABELS = {'Polarity 1', 'Polarity 2'};

CH_LIST        = 1:3;        % contacts to try; missing chN.mat files are skipped
NTRIG_EXPECTED = 10;         % pulses per amplitude level
PRETIME_MS     = 1.0;        % ms before stimulus in each epoch
POSTTIME_MS    = 5.0;        % ms after stimulus
ZOOM_XLIM      = [-0.5, 3.0];% ms, plot window
LINE_FREQ      = 60;         % Hz, passed to filter_periodic_artifact
LOW_AI_OVERRIDE = [];        % [] = lowest amplitude; or set an index (applies to all runs)

SAVE_FIGS      = false;      % save PNGs?
OUT_DIR        = fullfile(pwd, 'artifact_batch_out');
% ──────────────────────────────────────────────────────────────────────────────

stim_period_ms = 1 / 0.003096;                    % ~323 ms (3.096 Hz stim)
edge_thresh    = 0.01;

n_pol = numel(POLARITY_FOLDERS);
n_ch  = numel(CH_LIST);
results = cell(n_pol, max(CH_LIST));

if SAVE_FIGS && ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

for pi_ = 1:n_pol
    for ch = CH_LIST
        raw_path = fullfile(POLARITY_FOLDERS{pi_}, sprintf('ch%d.mat', ch));
        if ~isfile(raw_path)
            fprintf('[%s] ch%d: not found, skipping (%s)\n', POLARITY_LABELS{pi_}, ch, raw_path);
            continue;
        end
        fprintf('\n=== %s | ch%d ===\nLoading: %s\n', POLARITY_LABELS{pi_}, ch, raw_path);

        try
            res = process_contact(raw_path, NTRIG_EXPECTED, PRETIME_MS, POSTTIME_MS, ...
                                  stim_period_ms, edge_thresh, LINE_FREQ, LOW_AI_OVERRIDE);
        catch ME
            warning('[%s] ch%d failed: %s', POLARITY_LABELS{pi_}, ch, ME.message);
            continue;
        end
        res.polarity_label = POLARITY_LABELS{pi_};
        res.folder         = POLARITY_FOLDERS{pi_};
        res.ch             = ch;
        results{pi_, ch}   = res;

        fig = plot_contact(res, ZOOM_XLIM);
        if SAVE_FIGS
            [~, dsname] = fileparts(POLARITY_FOLDERS{pi_});
            exportgraphics(fig, fullfile(OUT_DIR, sprintf('pol%d_%s_ch%d.png', pi_, dsname, ch)), 'Resolution', 200);
        end
    end
end

%% Summary figure: mean before vs after, rows = polarity, cols = contact
fig_s = figure('Name', 'Summary: before vs after', 'NumberTitle', 'off', ...
               'Position', [50 50 350*n_ch 300*n_pol]);
tl = tiledlayout(n_pol, n_ch, 'TileSpacing', 'compact', 'Padding', 'compact');
for pi_ = 1:n_pol
    for ch = CH_LIST
        nexttile;
        if ch > size(results, 2) || isempty(results{pi_, ch})
            axis off; title(sprintf('%s ch%d (no data)', POLARITY_LABELS{pi_}, ch));
            continue;
        end
        r = results{pi_, ch};
        m = r.t_epoch >= ZOOM_XLIM(1) & r.t_epoch <= ZOOM_XLIM(2);
        hold on;
        plot(r.t_epoch(m), r.high_raw_avg(m),   'k-', 'LineWidth', 1.5);
        plot(r.t_epoch(m), r.high_clean_avg(m), '-', 'Color', [0.1 0.5 0.9], 'LineWidth', 1.5);
        xline(0, ':', 'Color', [0.5 0.5 0.5]);
        hold off; grid on; xlim(ZOOM_XLIM);
        title(sprintf('%s ch%d (%.3g mA)', POLARITY_LABELS{pi_}, ch, r.amp_nominal(r.high_ai)));
        if pi_ == n_pol, xlabel('Time (ms)'); end
        if ch == CH_LIST(1), ylabel('\muV'); end
        if pi_ == 1 && ch == CH_LIST(end), legend({'Before', 'After'}, 'Location', 'best'); end
    end
end
title(tl, 'Mean response before vs after template subtraction');
if SAVE_FIGS
    exportgraphics(fig_s, fullfile(OUT_DIR, 'summary.png'), 'Resolution', 200);
end

%% Save results
out_mat = fullfile(pwd, 'artifact_batch_results.mat');
save(out_mat, 'results', 'POLARITY_FOLDERS', 'POLARITY_LABELS', 'CH_LIST', '-v7.3');
fprintf('\nSaved results to %s\n', out_mat);

%% ── Local functions ──────────────────────────────────────────────────────────

function res = process_contact(raw_path, ntrig, pretime_ms, posttime_ms, ...
                               stim_period_ms, edge_thresh, line_freq, low_override)
    raw = load(raw_path);

    fs          = raw.fs;
    dt_ms       = 1000 / fs;
    amp_nominal = raw.amp_values_nominal(:)';
    n_amp       = numel(raw.ch_neural);
    fprintf('fs=%g Hz  n_amp=%d  amps: %s mA\n', fs, n_amp, num2str(amp_nominal, '%.3g '));

    min_gap_samp = round(0.9 * stim_period_ms / 1000 * fs);

    pre_samp  = round(pretime_ms  / 1000 * fs);
    post_samp = round(posttime_ms / 1000 * fs);
    n_epoch   = pre_samp + post_samp + 1;
    t_epoch   = (-pre_samp:post_samp) * dt_ms;   % ms, t=0 at stimulus

    epochs_all = cell(1, n_amp);
    for ai = 1:n_amp
        neural = raw.ch_neural{ai};
        [filtered_neural, ~] = filter_periodic_artifact(neural, fs, line_freq);
        trigs = find_triggers(raw.ch_Vinput{ai}, ntrig, min_gap_samp, edge_thresh);

        ep = zeros(n_epoch, ntrig);
        for ti = 1:ntrig
            s = trigs(ti) - pre_samp;
            e = trigs(ti) + post_samp;
            if s < 1 || e > numel(filtered_neural)
                ep(:, ti) = NaN;
            else
                ep(:, ti) = filtered_neural(s:e);
            end
        end
        epochs_all{ai} = ep * 1e6 / raw.gain_value;   % V -> uV, remove amp gain
    end

    % Template from lowest amplitude, subtract from highest
    [~, low_ai]  = min(amp_nominal);
    if ~isempty(low_override), low_ai = low_override; end
    [~, high_ai] = max(amp_nominal);

    low_epochs = epochs_all{low_ai};
    template   = mean(low_epochs, 2, 'omitnan');

    amp_scale = amp_nominal(high_ai) / amp_nominal(low_ai);
    fprintf('Template: %.3g mA (idx %d) | Subtract from: %.3g mA (idx %d) | scale %.4f\n', ...
        amp_nominal(low_ai), low_ai, amp_nominal(high_ai), high_ai, amp_scale);

    high_epochs = epochs_all{high_ai};
    high_clean  = high_epochs - amp_scale * template;   % implicit expansion over trials

    res = struct( ...
        't_epoch',        t_epoch, ...
        'amp_nominal',    amp_nominal, ...
        'low_ai',         low_ai, ...
        'high_ai',        high_ai, ...
        'amp_scale',      amp_scale, ...
        'epochs_all',     {epochs_all}, ...
        'low_epochs',     low_epochs, ...
        'template',       template, ...
        'high_epochs',    high_epochs, ...
        'high_clean',     high_clean, ...
        'high_raw_avg',   mean(high_epochs, 2, 'omitnan'), ...
        'high_clean_avg', mean(high_clean,  2, 'omitnan'));
end

function fig = plot_contact(r, zoom_xlim)
    m  = r.t_epoch >= zoom_xlim(1) & r.t_epoch <= zoom_xlim(2);
    nt = size(r.low_epochs, 2);
    gc = repmat(linspace(0.7, 0.1, nt)', 1, 3);
    col_template = [0.85 0.2 0.2];
    col_clean    = [0.1 0.5 0.9];

    fig = figure('Name', sprintf('%s ch%d', r.polarity_label, r.ch), 'NumberTitle', 'off', ...
                 'Position', [80 80 1100 750]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (1) low-amp epochs + template
    nexttile; hold on;
    for ti = 1:nt
        plot(r.t_epoch(m), r.low_epochs(m, ti), 'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
    end
    plot(r.t_epoch(m), r.template(m), 'Color', col_template, 'LineWidth', 2);
    xline(0, ':', 'Color', [0.5 0.5 0.5]); hold off; grid on; xlim(zoom_xlim);
    ylabel('\muV'); title(sprintf('Template epochs (%.3g mA)', r.amp_nominal(r.low_ai)));

    % (2) high-amp before
    nexttile; hold on;
    for ti = 1:nt
        plot(r.t_epoch(m), r.high_epochs(m, ti), 'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
    end
    plot(r.t_epoch(m), r.high_raw_avg(m), 'b-', 'LineWidth', 2);
    xline(0, ':', 'Color', [0.5 0.5 0.5]); hold off; grid on; xlim(zoom_xlim);
    title(sprintf('Before subtraction (%.3g mA)', r.amp_nominal(r.high_ai)));

    % (3) high-amp after
    nexttile; hold on;
    for ti = 1:nt
        plot(r.t_epoch(m), r.high_clean(m, ti), 'Color', [gc(ti,:) 0.5], 'LineWidth', 0.8);
    end
    plot(r.t_epoch(m), r.high_clean_avg(m), 'Color', col_clean, 'LineWidth', 2);
    xline(0, ':', 'Color', [0.5 0.5 0.5]); hold off; grid on; xlim(zoom_xlim);
    xlabel('Time (ms)'); ylabel('\muV'); title('After template subtraction');

    % (4) mean overlay
    nexttile; hold on;
    plot(r.t_epoch(m), r.high_raw_avg(m),   'k-', 'LineWidth', 1.8);
    plot(r.t_epoch(m), r.high_clean_avg(m), '-', 'Color', col_clean, 'LineWidth', 1.8);
    xline(0, ':', 'Color', [0.5 0.5 0.5]); hold off; grid on; xlim(zoom_xlim);
    xlabel('Time (ms)'); title('Mean before vs after'); legend({'Before', 'After'}, 'Location', 'best');

    title(tl, sprintf('%s  |  ch%d  |  %s', r.polarity_label, r.ch, r.folder), 'Interpreter', 'none');
end

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