% build_preproc_artifact_removed.m
%
% Processes all monophasic CAP datasets (all nerves with matching histology),
% removes the stim artifact using the low-amplitude template subtraction
% method from build_artifact_removed_df_SN/TIBIAL.m, and saves artifact-
% removed CAPData_ch*.mat files into a mirrored directory structure at:
%
%   /Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data_artifact_removed/
%
% The output mat files have identical fields to the originals, so
% cnap_analyzer.py can load them without modification.

%% ── Config ───────────────────────────────────────────────────────────────────
PREPROC_ROOT   = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data';
OUTPUT_ROOT    = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data_artifact_removed';
RAW_CACHE_ROOT = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data';

script_dir = fileparts(mfilename('fullpath'));
config_file = fullfile(script_dir, '..', 'local_config.yml');
if ~exist(config_file, 'file'); config_file = 'local_config.yml'; end

fid = fopen(config_file);
if fid == -1; error('Cannot open local_config.yml'); end
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

% Add 60 Hz denoising scripts
addpath(fullfile(script_dir, '..', 'external_dependencies', 'custom_denoising_scripts'));

%% ── Analysis parameters ──────────────────────────────────────────────────────
NTRIG_EXPECTED  = 10;
STIM_FREQ_HZ    = 3.096;
MIN_CV_MS       = 5.0;     % m/s — used for analysis window end

% Nerves and their expdates with matching histology
HIST_MAP = struct( ...
    'SN',     {{'20220531','20220720','20221108','20221201','20230207'}}, ...
    'tibial', {{'20220531','20220720','20221108','20221201','20230207'}}, ...
    'cVN',    {{'20220622','20221103','20221122','20230131'}} ...
);

%% ── Load metadata ────────────────────────────────────────────────────────────
meta = readtable(XLSX_PATH);
meta.expdate_str = arrayfun(@(x) num2str(x), meta.expdate_values, ...
                             'UniformOutput', false);

% Keep monophasic, valid polarity, matching histology
nerve_names = fieldnames(HIST_MAP);
keep = false(height(meta), 1);
for ni = 1:numel(nerve_names)
    nerve  = nerve_names{ni};
    dates  = HIST_MAP.(nerve);
    keep   = keep | (strcmp(meta.rec_nerve, nerve) & ...
                     strcmp(meta.phase_type, 'monophasic') & ...
                     ~isnan(meta.dataset_stim_polarity) & ...
                     ismember(meta.expdate_str, dates));
end
meta_proc = meta(keep, :);
fprintf('Processing %d datasets across all nerves\n', height(meta_proc));

%% ── Helpers ──────────────────────────────────────────────────────────────────
stim_period_ms  = 1000 / STIM_FREQ_HZ;
min_gap_samp_fn = @(fs) round(0.9 * stim_period_ms / 1000 * fs);
edge_thresh     = 0.01;
n_chs           = 3;

%% ── Main loop ────────────────────────────────────────────────────────────────
n_saved = 0;
n_skip  = 0;

for ri = 1:height(meta_proc)
    row       = meta_proc(ri, :);
    prefix    = row.file_prefix_outputs_values{1};
    stim_PW   = row.stim_PW;
    polarity  = row.dataset_stim_polarity;
    stim_dist = row.stim_rec_distance_center_to_center_mm;

    cache_dir   = fullfile(RAW_CACHE_ROOT, prefix);   % kept for existence check only
    preproc_dir = fullfile(PREPROC_ROOT,  prefix);
    out_dir     = fullfile(OUTPUT_ROOT,   prefix);

    if ~exist(preproc_dir, 'dir')
        fprintf('[SKIP %d/%d] %s — preproc dir missing\n', ...
            ri, height(meta_proc), prefix);
        n_skip = n_skip + 1;
        continue
    end

    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    fprintf('\n[%d/%d] %s  pol=%d  PW=%.3f ms  dist=%.1f mm\n', ...
        ri, height(meta_proc), prefix, polarity, stim_PW, stim_dist);

    for ch = 1:n_chs
        preproc_file = fullfile(preproc_dir, sprintf('CAPData_ch%d.mat', ch));
        out_file     = fullfile(out_dir,     sprintf('CAPData_ch%d.mat', ch));

        if ~exist(preproc_file, 'file')
            fprintf('  ch%d: missing preproc file, skipping\n', ch);
            continue
        end

        preproc = load(preproc_file);

        t_epoch       = preproc.t_window(:)';      % 1 x n_samp, ms
        n_samp_orig   = numel(t_epoch);
        amp_actual    = preproc.amp_values_actual(:)';   % 1 x n_amp
        n_amp_preproc = numel(amp_actual);

        %% Template = preprocessed CAP at the MINIMUM amplitude
        [~, low_ai]  = min(amp_actual);
        [~, high_ai] = max(amp_actual);
        fprintf('  ch%d: template=%.3g mA (col %d)  max=%.3g mA (col %d)\n', ...
            ch, amp_actual(low_ai), low_ai, amp_actual(high_ai), high_ai);

        low_cap  = preproc.ch_neural_avg(:, low_ai);   % n_samp x 1
        high_cap = preproc.ch_neural_avg(:, high_ai);  % for diagnostics

        %% Subtract scaled template from every amplitude column
        ch_neural_avg_clean = zeros(n_samp_orig, n_amp_preproc);
        for ai_p = 1:n_amp_preproc
            amp_scale = amp_actual(ai_p) / amp_actual(low_ai);
            ch_neural_avg_clean(:, ai_p) = preproc.ch_neural_avg(:, ai_p) ...
                                            - amp_scale * low_cap;
        end

        %% Diagnostic plot for ch1 only (COMMENTED OUT)
        if ch == 1
            zoom = t_epoch >= -0.5 & t_epoch <= 3.0;
            amp_scale_high = amp_actual(high_ai) / amp_actual(low_ai);
            high_clean = ch_neural_avg_clean(:, high_ai);

            figure('Name', sprintf('Diag: %s ch1', prefix), ...
                   'NumberTitle', 'off', 'Position', [50 50 1100 500]);

            subplot(1,3,1);
            plot(t_epoch(zoom), low_cap(zoom), 'r-', 'LineWidth', 2);
            xline(0,'k:');
            title(sprintf('Template  %.3g mA', amp_actual(low_ai)), 'FontSize', 11);
            xlabel('Time (ms)'); ylabel('µV'); grid on;

            subplot(1,3,2);
            plot(t_epoch(zoom), high_cap(zoom), 'k--', 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Before  %.3g mA', amp_actual(high_ai)));
            hold on;
            plot(t_epoch(zoom), high_clean(zoom), 'b-', 'LineWidth', 2, ...
                 'DisplayName', 'After subtraction');
            xline(0,'k:'); legend('FontSize', 9);
            title('Before vs After (max amp)', 'FontSize', 11);
            xlabel('Time (ms)'); ylabel('µV'); grid on;

            subplot(1,3,3);
            cmap = parula(n_amp_preproc);
            hold on;
            for ai_p = 1:n_amp_preproc
                plot(t_epoch(zoom), ch_neural_avg_clean(zoom, ai_p), ...
                     'Color', cmap(ai_p,:), 'LineWidth', 1.2);
            end
            xline(0,'k:');
            title('All amplitudes — artifact removed', 'FontSize', 11);
            xlabel('Time (ms)'); ylabel('µV'); grid on;
            cb = colorbar; cb.Label.String = 'Amp (mA)';
            colormap(parula); clim([amp_actual(1) amp_actual(end)]);

            sgtitle(sprintf('%s  dist=%.1fmm  pol=%d  PW=%.3fms', ...
                prefix, stim_dist, polarity, stim_PW), 'FontSize', 12);
            drawnow;
        end

        %% Save: copy original mat, replace ch_neural_avg with artifact-removed signal
        out_data = preproc;
        out_data.ch_neural_avg = ch_neural_avg_clean;

        % Recompute AUC for artifact-removed signal
        art_end_idx = find(t_epoch >= 0, 1, 'first');
        ana_end_idx = find(t_epoch >= stim_dist / MIN_CV_MS, 1, 'first');
        AUC_clean = zeros(1, n_amp_preproc);
        if ~isempty(art_end_idx) && ~isempty(ana_end_idx) && ana_end_idx > art_end_idx
            for ai_p = 1:n_amp_preproc
                AUC_clean(ai_p) = trapz( ...
                    t_epoch(art_end_idx:ana_end_idx), ...
                    abs(ch_neural_avg_clean(art_end_idx:ana_end_idx, ai_p)));
            end
        end
        out_data.AUC = AUC_clean;

        % Save to output directory
        save(out_file, '-struct', 'out_data');
        fprintf('  ch%d: saved → %s\n', ch, out_file);
        n_saved = n_saved + 1;
    end
end

fprintf('\nDone. Saved %d files, skipped %d datasets.\n', n_saved, n_skip);
fprintf('Output root: %s\n', OUTPUT_ROOT);

%% ── Local functions ──────────────────────────────────────────────────────────
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
