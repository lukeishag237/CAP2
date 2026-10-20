% compare_amplitudes.m
%
% Loads processed CAP data from two folders and creates three visualization plots:
% 1. All amplitudes overlaid from Folder 1
% 2. All amplitudes overlaid from Folder 2
% 3. Combined overlay comparing both folders at each amplitude level
%
% Just edit the user settings below and press Run!

% ── USER SETTINGS ─────────────────────────────────────────────────────────────
FOLDER1 = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220525_7';
FOLDER2 = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data/20220525_8';
CH_NUM  = 2;
% ──────────────────────────────────────────────────────────────────────────────

compare_amplitudes_func(FOLDER1, FOLDER2, CH_NUM);

function compare_amplitudes_func(folder1, folder2, ch_num)
    
    ZOOM_XLIM      = [-0.5, 3.0];   % ms
    PRETIME_MS     = 1.0;
    POSTTIME_MS    = 5.0;
    
    fprintf('Loading data from folder 1: %s\n', folder1);
    [amp1, epochs1, t1, fs1] = load_amplitude_data(folder1, ch_num, PRETIME_MS, POSTTIME_MS);
    n_amp1 = length(amp1);
    
    fprintf('Loading data from folder 2: %s\n', folder2);
    [amp2, epochs2, t2, fs2] = load_amplitude_data(folder2, ch_num, PRETIME_MS, POSTTIME_MS);
    n_amp2 = length(amp2);
    
    zoom_mask1 = t1 >= ZOOM_XLIM(1) & t1 <= ZOOM_XLIM(2);
    zoom_mask2 = t2 >= ZOOM_XLIM(1) & t2 <= ZOOM_XLIM(2);
    
    cmap1 = viridis(n_amp1);
    cmap2 = plasma(n_amp2);
    
    % ========================================================================
    % Plot 1: All amplitudes from Folder 1
    % ========================================================================
    fig1 = figure('Name', sprintf('Folder 1: All amplitudes — Ch%d', ch_num), ...
                  'NumberTitle', 'off', 'Position', [50 80 1200 600]);
    ax1 = axes('Parent', fig1);
    hold(ax1, 'on');
    
    for ai = 1:n_amp1
        plot(ax1, t1(zoom_mask1), epochs1(zoom_mask1, ai), ...
            'Color', cmap1(ai, :), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%.3g mA', amp1(ai)));
    end
    
    xline(ax1, 0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    hold(ax1, 'off');
    
    ax1.XLabel.String = 'Time (ms)';
    ax1.YLabel.String = 'Amplitude (µV)';
    ax1.Title.String = sprintf('Folder 1: All amplitudes overlaid (Ch%d)', ch_num);
    ax1.Title.FontSize = 13;
    ax1.XLabel.FontSize = 12;
    ax1.YLabel.FontSize = 12;
    ax1.XLim = ZOOM_XLIM;
    
    legend(ax1, 'Location', 'eastoutside', 'FontSize', 8, 'NumColumns', 2);
    
    file1 = fullfile(folder1, sprintf('compare_amplitudes_folder1_ch%d.png', ch_num));
    exportgraphics(fig1, file1, 'Resolution', 150);
    fprintf('Saved: %s\n', file1);
    
    % ========================================================================
    % Plot 2: All amplitudes from Folder 2
    % ========================================================================
    fig2 = figure('Name', sprintf('Folder 2: All amplitudes — Ch%d', ch_num), ...
                  'NumberTitle', 'off', 'Position', [100 80 1200 600]);
    ax2 = axes('Parent', fig2);
    hold(ax2, 'on');
    
    for ai = 1:n_amp2
        plot(ax2, t2(zoom_mask2), epochs2(zoom_mask2, ai), ...
            'Color', cmap2(ai, :), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%.3g mA', amp2(ai)));
    end
    
    xline(ax2, 0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    hold(ax2, 'off');
    
    ax2.XLabel.String = 'Time (ms)';
    ax2.YLabel.String = 'Amplitude (µV)';
    ax2.Title.String = sprintf('Folder 2: All amplitudes overlaid (Ch%d)', ch_num);
    ax2.Title.FontSize = 13;
    ax2.XLabel.FontSize = 12;
    ax2.YLabel.FontSize = 12;
    ax2.XLim = ZOOM_XLIM;
    
    legend(ax2, 'Location', 'eastoutside', 'FontSize', 8, 'NumColumns', 2);
    
    file2 = fullfile(folder2, sprintf('compare_amplitudes_folder2_ch%d.png', ch_num));
    exportgraphics(fig2, file2, 'Resolution', 150);
    fprintf('Saved: %s\n', file2);
    
    % ========================================================================
    % Plot 3: Combined comparison at each amplitude level
    % ========================================================================
    % Get all unique amplitude values
    all_amps = unique([amp1; amp2]);
    all_amps = sort(all_amps);
    n_subplots = length(all_amps);
    
    % Grid layout
    ncols = min(n_subplots, 5);
    nrows = ceil(n_subplots / ncols);
    
    % Calculate global y limits
    y_all = [epochs1(zoom_mask1, :); epochs2(zoom_mask2, :)];
    yrange = max(y_all(:)) - min(y_all(:));
    ylims = [min(y_all(:)) - 0.05*yrange, max(y_all(:)) + 0.05*yrange];
    
    fig3 = figure('Name', sprintf('Amplitude Comparison — Ch%d', ch_num), ...
                  'NumberTitle', 'off', 'Position', [150 100 400*ncols 300*nrows]);
    
    ax_all = gobjects(nrows, ncols);
    
    for idx = 1:n_subplots
        amp_val = all_amps(idx);
        r = ceil(idx / ncols);
        c = mod(idx - 1, ncols) + 1;
        
        ax = subplot(nrows, ncols, idx);
        ax_all(r, c) = ax;
        hold(ax, 'on');
        
        % Plot from folder 1 if this amplitude exists
        [~, ai1] = ismember(amp_val, amp1);
        if ai1 > 0
            plot(ax, t1(zoom_mask1), epochs1(zoom_mask1, ai1), ...
                'Color', '#1F77B4', 'LineWidth', 1.5, 'DisplayName', 'Folder 1');
        end
        
        % Plot from folder 2 if this amplitude exists
        [~, ai2] = ismember(amp_val, amp2);
        if ai2 > 0
            plot(ax, t2(zoom_mask2), epochs2(zoom_mask2, ai2), ...
                'Color', '#FF7F0E', 'LineWidth', 1.5, 'DisplayName', 'Folder 2');
        end
        
        xline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 1, ...
            'HandleVisibility', 'off');
        hold(ax, 'off');
        
        ax.Title.String = sprintf('%.3g mA', amp_val);
        ax.Title.FontSize = 10;
        ax.Title.FontWeight = 'bold';
        ax.XLim = ZOOM_XLIM;
        ax.YLim = ylims;
        
        if c == 1
            ax.YLabel.String = 'µV';
            ax.YLabel.FontSize = 9;
        end
        
        if r == nrows
            ax.XLabel.String = 'Time (ms)';
            ax.XLabel.FontSize = 9;
        end
        
        if idx == 1
            legend(ax, 'Location', 'northeast', 'FontSize', 8);
        end
    end
    
    % Add overall title
    sgtitle(fig3, sprintf('Amplitude Comparison: Folder 1 vs Folder 2 (Ch%d)', ch_num), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Save comparison plot
    [parent_dir, ~] = fileparts(folder1);
    file3 = fullfile(parent_dir, sprintf('compare_amplitudes_overlay_ch%d.png', ch_num));
    exportgraphics(fig3, file3, 'Resolution', 150);
    fprintf('Saved: %s\n', file3);
    
    fprintf('Done!\n');
    
end

% ============================================================================
% LOCAL FUNCTIONS
% ============================================================================

function [amp_values, avg_epochs, t_epoch, fs] = load_amplitude_data(folder_path, ch_num, pretime_ms, posttime_ms)
    
    % Try to load from MAT file
    raw_path = fullfile(folder_path, sprintf('ch%d.mat', ch_num));
    
    if ~isfile(raw_path)
        error('No data file found: %s', raw_path);
    end
    
    fprintf('Loading: %s\n', raw_path);
    raw = load(raw_path);
    
    fs = raw.fs;
    dt_ms = 1000 / fs;
    amp_values = raw.amp_values_nominal(:);
    n_amp = length(raw.ch_neural);
    
    fprintf('n_amp=%d  amps: ', n_amp);
    fprintf('%.3g ', amp_values);
    fprintf('mA\n');
    
    % Process data
    stim_period_ms = 1000 / 3.096;
    min_gap = round(0.9 * stim_period_ms / 1000 * fs);
    edge_thresh = 0.01;
    
    pre_samp = round(pretime_ms / dt_ms);
    post_samp = round(posttime_ms / dt_ms);
    n_epoch = pre_samp + post_samp + 1;
    t_epoch = (-pre_samp:post_samp)' * dt_ms;
    
    ntrig_expected = 10;
    avg_epochs = zeros(n_epoch, n_amp);
    
    for ai = 1:n_amp
        neural = raw.ch_neural{ai};
        vinput = raw.ch_Vinput{ai};
        
        trigs = find_triggers(vinput, ntrig_expected, min_gap, edge_thresh);
        
        ep = zeros(n_epoch, ntrig_expected);
        for ti = 1:length(trigs)
            s = trigs(ti) - pre_samp;
            e = trigs(ti) + post_samp;
            if s >= 1 && e <= length(neural)
                trial = neural(s:e) * 1e6 / raw.gain_value;
                trial = trial - mean(trial(t_epoch < 0));  % baseline correct
                ep(:, ti) = trial;
            end
        end
        avg_epochs(:, ai) = mean(ep, 2);
    end
    
end

function trig = find_triggers(vinput_col, ntrig, min_gap, edge_thr)
    
    abs_d = diff(abs(vinput_col));
    trig = zeros(1, ntrig);
    
    % Find first trigger
    idx = find(abs_d > edge_thr, 1, 'first');
    if isempty(idx)
        trig = [];
        return;
    end
    trig(1) = idx;
    
    % Find remaining triggers
    for p = 2:ntrig
        s = trig(p-1) + min_gap;
        if s >= length(abs_d)
            break;
        end
        r = find(abs_d(s:end) > edge_thr, 1, 'first');
        if ~isempty(r)
            trig(p) = s - 1 + r;
        else
            break;
        end
    end
    
end

function cmap = viridis(n)
    % Generate viridis colormap with n colors
    if n == 1
        cmap = [0.267004, 0.004874, 0.329415];
    else
        x = linspace(0, 1, n)';
        % Viridis colormap approximation
        cmap = interp1([0 0.25 0.5 0.75 1], ...
            [0.267004 0.282623 0.253935 0.745693 0.993248; ...
             0.004874 0.140461 0.528596 0.905164 0.906157; ...
             0.329415 0.615762 0.715467 0.143080 0.143640]', x);
    end
end

function cmap = plasma(n)
    % Generate plasma colormap with n colors
    if n == 1
        cmap = [0.050383, 0.029803, 0.527975];
    else
        x = linspace(0, 1, n)';
        % Plasma colormap approximation
        cmap = interp1([0 0.25 0.5 0.75 1], ...
            [0.050383 0.282623 0.798340 0.940015 0.940015; ...
             0.029803 0.140461 0.009127 0.975158 0.975158; ...
             0.527975 0.615762 0.914313 0.131326 0.131326]', x);
    end
end
