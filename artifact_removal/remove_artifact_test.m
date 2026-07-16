baseFolder = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data';

% Get all directories matching the pattern 20220720_*
dirs = dir(fullfile(baseFolder, '20220720_*'));

% Keep only actual directories (exclude files)
dirs = dirs([dirs.isdir]);

fprintf('Found %d directories\n', numel(dirs));

for d = 1:numel(dirs)
    thisDir = fullfile(dirs(d).folder, dirs(d).name);
    fprintf('\nEntering directory: %s\n', thisDir);

    % Find all .mat files inside this directory
    matFiles = dir(fullfile(thisDir, '*.mat'));
    fprintf('  Found %d .mat files\n', numel(matFiles));

    for k = 1:numel(matFiles)
        fname = fullfile(matFiles(k).folder, matFiles(k).name);
        fprintf('    Loading %s\n', fname);

        S = load(fname);

        % ---- Extract waveform and time vector ----
        x = S.ch_neural_avg;      % 32301 x 11
        disp(size(x))
        t = S.t_window;           % 1 x 32301, in ms
        
        % ---- Define time window ----
        t_start = 0;      % ms
        t_end   = 10;     % ms
        
        % ---- Find indices for the window ----
        idx = (t >= t_start) & (t <= t_end);
        
        % ---- Plot each column (each trial) ----
        figure;
        plot(t(idx), x(idx, :), 'LineWidth', 1.2);
        title(sprintf('%s / %s', dirs(d).name, matFiles(k).name), ...
              'Interpreter', 'none');
        xlabel('Time (ms)');
        ylabel('Amplitude (a.u.)');
        grid on;

    end
end
