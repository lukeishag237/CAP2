% build_artifact_removed_df.m
%
% For all SN monophasic datasets: load raw cache, subtract stimulus artifact
% (template = mean of 10 trials at lowest amplitude), re-extract peak
% amplitude, peak time, and AUC, then save a table as a .mat file that
% mirrors stored_df.pkl.
%
% Output: <local_output_path>/stored_df_artifact_removed_SN.mat
%         Contains a struct array 'records' with one entry per channel per
%         dataset, plus a table 'T' with the same columns as stored_df.pkl.

%% ── Config ──────────────────────────────────────────────────────────────────
local_config_file = 'local_config.yml';
fid = fopen(local_config_file);
raw_data_path = ''; preprocessed_path = ''; local_output_path = '';
line_i = fgetl(fid);
while ischar(line_i)
    if contains(line_i,'raw_data_path: ')
        raw_data_path = strtrim(line_i(length('raw_data_path: ')+2:end-1));
    elseif contains(line_i,'preprocessed_CAP_path: ')
        preprocessed_path = strtrim(line_i(length('preprocessed_CAP_path: ')+2:end-1));
    elseif contains(line_i,'local_output_path: ')
        local_output_path = strtrim(line_i(length('local_output_path: ')+2:end-1));
    end
    line_i = fgetl(fid);
end
fclose(fid);

RAW_CACHE_ROOT = fullfile(raw_data_path, '..', 'mat_cache_of_raw_data');
XLSX_PATH      = fullfile(raw_data_path, 'AnalysisScript_CAPs_Batch_parameters.xlsx');
OUTPUT_PATH    = fullfile(local_output_path, 'stored_df_artifact_removed_SN.mat');

%% ── Analysis parameters ─────────────────────────────────────────────────────
NTRIG_EXPECTED = 10;
MIN_CV_MS      = 5.0;    % m/s  (used for analysis window end)
STIM_FREQ_HZ   = 3.096;  % Hz
TEMPLATE_AMP_MA = 0.08;   % mA — use amplitude closest to this for the artifact template

% Channel offsets for 3-contact cuff [mm]
CH_OFFSET = [-1.5, 0.0, 1.5];

%% ── Load metadata ────────────────────────────────────────────────────────────
meta = readtable(XLSX_PATH);
meta.expdate_str = arrayfun(@(x) num2str(x), meta.expdate_values, ...
                             'UniformOutput', false);

% Keep only SN monophasic with valid polarity AND matching histology
% Histology mapping (from compare_cnap_and_histology.py / 3d_view.py):
%   SN expdates with histology: 20220531, 20220720, 20221108, 20221201, 20230207
SN_HIST_EXPDATES = {'20220531','20220720','20221108','20221201','20230207'};

keep = strcmp(meta.rec_nerve, 'SN') & ...
       strcmp(meta.phase_type, 'monophasic') & ...
       ~isnan(meta.dataset_stim_polarity) & ...
       ismember(meta.expdate_str, SN_HIST_EXPDATES);
meta_sn = meta(keep, :);

fprintf('Processing %d SN monophasic datasets\n', height(meta_sn));

%% ── Helper parameters ────────────────────────────────────────────────────────
stim_period_ms  = 1000 / STIM_FREQ_HZ;
min_gap_samp_fn = @(fs) round(0.9 * stim_period_ms / 1000 * fs);
edge_thresh     = 0.01;
n_chs           = 3;

%% ── Main loop ────────────────────────────────────────────────────────────────
all_rows = {};   % will become table rows

for ri = 1:height(meta_sn)
    row        = meta_sn(ri, :);
    prefix     = row.file_prefix_outputs_values{1};
    expdate    = row.expdate_str{1};
    stim_PW    = row.stim_PW;
    polarity   = row.dataset_stim_polarity;
    stim_dist  = row.stim_rec_distance_center_to_center_mm;
    fiber_type = row.target_fiber_types{1};

    cache_dir = fullfile(RAW_CACHE_ROOT, prefix);
    if ~exist(cache_dir, 'dir')
        fprintf('  [SKIP] cache not found: %s\n', cache_dir);
        continue
    end

    fprintf('\n[%d/%d] %s  pol=%d  PW=%.3f ms  dist=%.1f mm\n', ...
        ri, height(meta_sn), prefix, polarity, stim_PW, stim_dist);

    % Process each channel
    for ch = 1:n_chs
        cache_file = fullfile(cache_dir, sprintf('ch%d.mat', ch));
        if ~exist(cache_file, 'file')
            fprintf('  ch%d: cache file missing, skipping\n', ch);
            continue
        end

        raw = load(cache_file);
        fs  = raw.fs;
        dt  = 1/fs;   % seconds
        dt_ms = dt * 1000;
        n_amp = numel(raw.ch_neural);
        amp_nominal = raw.amp_values_nominal(:)';
        min_gap = min_gap_samp_fn(fs);

        % ── Slice all epochs ────────────────────────────────────────────────
        % Epoch window: 2 ms pre, enough post to cover slowest fiber
        pretime_ms  = 2.0;
        posttime_ms = stim_dist / 0.5 + 2.0;  % cover slowest fiber + buffer
        pre_samp    = round(pretime_ms  / dt_ms);
        post_samp   = round(posttime_ms / dt_ms);
        n_epoch     = pre_samp + post_samp + 1;
        t_epoch     = (-pre_samp:post_samp) * dt_ms;   % ms, t=0 at stim

        epochs_all = zeros(n_epoch, NTRIG_EXPECTED, n_amp);

        for ai = 1:n_amp
            neural = raw.ch_neural{ai};
            vinput = raw.ch_Vinput{ai};
            trigs  = find_triggers(vinput, NTRIG_EXPECTED, min_gap, edge_thresh);

            for ti = 1:NTRIG_EXPECTED
                s = trigs(ti) - pre_samp;
                e = trigs(ti) + post_samp;
                if s >= 1 && e <= numel(neural)
                    ep = neural(s:e) * 1e6 / raw.gain_value;  % µV
                    % Pre-stim baseline correction
                    pre_mask_ep = t_epoch < 0;
                    ep = ep - mean(ep(pre_mask_ep));
                    epochs_all(:, ti, ai) = ep;
                end
            end
        end

        % ── Artifact template: amplitude closest to TEMPLATE_AMP_MA ─────────
        [~, low_ai]  = min(abs(amp_nominal - TEMPLATE_AMP_MA));
        [~, high_ai] = max(amp_nominal);
        fprintf('  ch%d: template from %.3g mA (index %d), subtract up to %.3g mA (index %d)\n', ...
            ch, amp_nominal(low_ai), low_ai, amp_nominal(high_ai), high_ai);

        low_avg  = mean(epochs_all(:, :, low_ai), 2);   % n_epoch x 1

        % ── Subtract template from every amplitude (scale by amp ratio) ─────
        epochs_clean = zeros(n_epoch, NTRIG_EXPECTED, n_amp);
        for ai = 1:n_amp
            amp_scale = amp_nominal(ai) / amp_nominal(low_ai);
            for ti = 1:NTRIG_EXPECTED
                epochs_clean(:, ti, ai) = epochs_all(:, ti, ai) - amp_scale * low_avg;
            end
        end

        % ── Average clean epochs (median) ────────────────────────────────────
        ch_neural_avg_clean = median(epochs_clean, 2);   % n_epoch x 1 x n_amp
        ch_neural_avg_clean = squeeze(ch_neural_avg_clean);  % n_epoch x n_amp

        % Flip polarity 2 back to original orientation
        if polarity == 2
            ch_neural_avg_clean = -ch_neural_avg_clean;
        end

        % ── Plot: Fig 1 (template) + Fig 3 (before vs after mean) ────────────
        if ch == 1
        zoom_mask   = t_epoch >= -0.5 & t_epoch <= 3.0;
        col_template = [0.85 0.2 0.2];
        col_clean    = [0.1  0.5 0.9];
        gc           = repmat(linspace(0.7, 0.2, NTRIG_EXPECTED)', 1, 3);
        fig_title_base = sprintf('%s  ch%d  dist=%.1fmm  pol%d', ...
            prefix, ch, stim_dist, polarity);

        % Figure 1: low-amp epochs + template
        low_epochs = squeeze(epochs_all(:, :, low_ai));
        if polarity == 2; low_epochs = -low_epochs; end
        low_avg_disp = mean(low_epochs, 2);

        figure('Name', sprintf('Template — %s ch%d', prefix, ch), ...
               'NumberTitle', 'off', 'Position', [60 120 750 420]);
        hold on;
        for ti = 1:NTRIG_EXPECTED
            plot(t_epoch(zoom_mask), low_epochs(zoom_mask, ti), ...
                 'Color', [gc(ti,:) 0.45], 'LineWidth', 0.8);
        end
        plot(t_epoch(zoom_mask), low_avg_disp(zoom_mask), ...
             'Color', col_template, 'LineWidth', 2.0);
        xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        hold off;
        xlabel('Time (ms)', 'FontSize', 11);
        ylabel('Amplitude (µV)', 'FontSize', 11);
        title(sprintf('Artifact template  %.3g mA — %s', ...
            amp_nominal(low_ai), fig_title_base), 'FontSize', 11);
        legend({'Trials','Mean template'}, 'FontSize', 9, 'Location', 'northeast');
        xlim([-0.5 3.0]); grid on;

        % Figure 3: mean before vs mean after at highest amplitude
        high_raw_avg   = mean(squeeze(epochs_all(:, :, high_ai)), 2);
        high_clean_avg = ch_neural_avg_clean(:, high_ai);
        % Flip raw for display if polarity 2
        if polarity == 2; high_raw_avg = -high_raw_avg; end

        figure('Name', sprintf('Before vs After — %s ch%d', prefix, ch), ...
               'NumberTitle', 'off', 'Position', [80 140 750 420]);
        hold on;
        plot(t_epoch(zoom_mask), high_raw_avg(zoom_mask),   'k-',  'LineWidth', 1.8, ...
             'DisplayName', sprintf('Before  %.3g mA', amp_nominal(high_ai)));
        plot(t_epoch(zoom_mask), high_clean_avg(zoom_mask), ...
             'Color', col_clean, 'LineWidth', 1.8, 'DisplayName', 'After subtraction');
        xline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        hold off;
        xlabel('Time (ms)', 'FontSize', 11);
        ylabel('Amplitude (µV)', 'FontSize', 11);
        title(sprintf('Before vs After — %s', fig_title_base), 'FontSize', 11);
        legend('FontSize', 10);
        xlim([-0.5 3.0]); grid on;

        drawnow;
        end % if ch == 1

        % ── Select max amplitude column ───────────────────────────────────────
        signal   = ch_neural_avg_clean(:, high_ai);
        amp_used = amp_nominal(high_ai);

        % ── Find peak in analysis window [artifact_end, dist/min_cv] ─────────
        artifact_end_ms = stim_PW + 0.17;
        analysis_end_ms = stim_dist / MIN_CV_MS;

        art_end_idx = find(t_epoch >= artifact_end_ms, 1, 'first');
        ana_end_idx = find(t_epoch >= analysis_end_ms, 1, 'first');
        if isempty(art_end_idx) || isempty(ana_end_idx) || ana_end_idx <= art_end_idx
            fprintf('  ch%d: window invalid, skipping\n', ch);
            continue
        end

        [peak_amp, peak_rel] = max(signal(art_end_idx:ana_end_idx));
        peak_idx  = art_end_idx + peak_rel - 1;
        peak_time = t_epoch(peak_idx);

        % AUC (trapz of abs signal in analysis window)
        auc = trapz(t_epoch(art_end_idx:ana_end_idx), ...
                    abs(signal(art_end_idx:ana_end_idx)));

        % ── Corrected stim-rec distance (per-channel offset) ─────────────────
        dist_corrected = stim_dist + CH_OFFSET(ch);

        % ── Store row ─────────────────────────────────────────────────────────
        all_rows{end+1} = struct( ...
            'expdate',              expdate, ...
            'rec_nerve',            'SN', ...
            'target_fiber_types',   fiber_type, ...
            'ch_num_neural',        ch, ...
            'stim_rec_dist',        dist_corrected, ...
            'peak_amplitude',       peak_amp, ...
            'peak_time',            peak_time, ...
            'auc',                  auc, ...
            'dataset_stim_polarity',polarity, ...
            't_window',             t_epoch, ...
            'ch_neural_avg',        signal' ...   % row vector for table compat
        );

        fprintf('  ch%d: peak=%.1f µV  t=%.3f ms  auc=%.2f\n', ...
            ch, peak_amp, peak_time, auc);
    end
end

%% ── Build table and save ─────────────────────────────────────────────────────
if isempty(all_rows)
    error('No rows collected — check cache paths and metadata');
end

% Convert struct array to table
n = numel(all_rows);
T = table();
T.expdate               = cellfun(@(r) r.expdate,              all_rows, 'UniformOutput', false)';
T.rec_nerve             = repmat({'SN'}, n, 1);
T.target_fiber_types    = cellfun(@(r) r.target_fiber_types,   all_rows, 'UniformOutput', false)';
T.ch_num_neural         = cellfun(@(r) r.ch_num_neural,        all_rows)';
T.stim_rec_dist         = cellfun(@(r) r.stim_rec_dist,        all_rows)';
T.peak_amplitude        = cellfun(@(r) r.peak_amplitude,       all_rows)';
T.peak_time             = cellfun(@(r) r.peak_time,            all_rows)';
T.auc                   = cellfun(@(r) r.auc,                  all_rows)';
T.dataset_stim_polarity = cellfun(@(r) r.dataset_stim_polarity,all_rows)';

fprintf('\nCollected %d rows\n', n);
disp(T(:, {'expdate','ch_num_neural','stim_rec_dist','peak_amplitude','peak_time','auc'}));

save(OUTPUT_PATH, 'T', 'all_rows');
fprintf('Saved to: %s\n', OUTPUT_PATH);

%% ── Local functions ──────────────────────────────────────────────────────────
function trig = find_triggers(vinput_col, ntrig, min_gap, edge_thr)
    abs_d  = diff(abs(vinput_col));
    trig   = zeros(1, ntrig);
    trig(1)= find(abs_d > edge_thr, 1, 'first');
    for p  = 2:ntrig
        s  = trig(p-1) + min_gap;
        r  = find(abs_d(s:end) > edge_thr, 1, 'first');
        trig(p) = s - 1 + r;
    end
end
