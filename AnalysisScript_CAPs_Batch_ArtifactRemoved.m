% AnalysisScript_CAPs_Batch_ArtifactRemoved.m
%
% Batch script to process datasets with artifact removal.
% Uses AnalysisScript_CAPs_ArtifactRemoved instead of AnalysisScript_CAPs_Wrapper.

clear
close all
format compact

scrnsz = get(0,'screensize');

% Add paths
addpath('external_dependencies/adinstruments_sdk_matlab-master');
addpath('external_dependencies/custom_denoising_scripts');

%% Get paths from YAML config
local_config_file = 'local_config.yml';
fid = fopen(local_config_file);
line_i = fgetl(fid);
raw_data_path = '';
preprocessed_CAP_path = '';
local_output_path = '';
while ischar(line_i)
    if contains(line_i,'raw_data_path: ')
        raw_data_path = line_i(length('raw_data_path: ')+2:end-1);
    elseif contains(line_i,'preprocessed_CAP_path: ')
        preprocessed_CAP_path = line_i(length('preprocessed_CAP_path: ')+2:end-1);
    elseif contains(line_i,'local_output_path: ')
        local_output_path = line_i(length('local_output_path: ')+2:end-1);
    end
    line_i = fgetl(fid);
end
fclose(fid);
assert(~strcmp(raw_data_path,''), 'raw_data_path not found in local_config.yml');
clearvars -except local_config_file raw_data_path preprocessed_CAP_path local_output_path

%% Artifact removed output path
artifact_removed_output_path = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/preprocessed_CAP_data_artifact_removed';

%% Skip if already exists flag
skip_if_exists = true;  % Set to false to force reprocessing of all datasets

%% Define analysis parameters
CAP_avg_type            = 2;                    % 1 for mean; 2 for median
subtractnoise           = 1;                    % Subtract noise using AUC for t < 0
zerobaseline            = 1;                    % Zero baseline
digitalLPF              = 0;                    % Digital LPF (not used in artifact removal version)
digitalLPF_fc_bandpass  = 0;
digitaloffsetremoval    = 1;

%% Load metadata
xlsx_path = fullfile(raw_data_path, 'AnalysisScript_CAPs_Batch_parameters.xlsx');
meta = readtable(xlsx_path);
meta.expdate_str = arrayfun(@(x) num2str(x), meta.expdate_values, ...
                             'UniformOutput', false);

% Keep only monophasic with valid polarity, OR cVN datasets
% cVN expdates: 20220622, 20221103, 20221122, 20230131
cVN_expdates = [20220622, 20221103, 20221122, 20230131];
keep = (strcmp(meta.phase_type, 'monophasic') & ~isnan(meta.dataset_stim_polarity)) | ...
       ismember(meta.expdate_values, cVN_expdates);
meta_proc = meta(keep, :);

fprintf('Processing %d datasets with artifact removal\n', height(meta_proc));

%% Process each dataset
for ri = 1:height(meta_proc)
    row = meta_proc(ri, :);
    
    expdate = row.expdate_values;
    file_prefix_metadata = row.file_prefix_metadata{1};
    file_prefix_adicht = row.file_prefix_adicht{1};
    
    % Use artifact-removed output path
    file_prefix_outputs = fullfile(artifact_removed_output_path, row.file_prefix_outputs_values{1});
    
    % Full paths
    file_prefix_metadata_full = fullfile(raw_data_path, file_prefix_metadata);
    file_prefix_adicht_full = fullfile(file_prefix_metadata_full, file_prefix_adicht);
    
    fprintf('\n[%d/%d] Processing %s (expdate=%s)\n', ...
        ri, height(meta_proc), file_prefix_outputs, num2str(expdate));
    
    % Load trial info from TestPulseData
    trialnum_start = row.trialnum_start_values;
    % trialnum_end is not in the xlsx - we use the trial_values from the TestPulseData file
    
    % Determine channel settings based on electrode type and monopolar recording
    ch_num_Vinput = 4;
    if (row.dataset_electrodetype == 1) || (row.dataset_monopolarrec == 0)
        ch_num_neural_values = 1;
    elseif (row.dataset_electrodetype == 2) && (row.dataset_monopolarrec == 1)
        ch_num_neural_values = [1 2];
    elseif (row.dataset_electrodetype == 3) && (row.dataset_monopolarrec == 1)
        ch_num_neural_values = [1 2 3];
    else
        error('Unexpected electrode type and monopolar recording combination');
    end
    
    pulse_width_us = row.stim_PW * 1000;
    amp_file = sprintf('%s/TestPulseData_Trial%d_ch%d.mat', ...
        file_prefix_metadata_full, trialnum_start, ch_num_neural_values(1));
    
    if exist(amp_file, 'file')
        load(amp_file, 'amp_values_nominal', 'trial_values');
        
        % Set actual amplitudes (no attenuation correction)
        amp_values_actual = amp_values_nominal;
        
        % Set up stimulation distance and other parameters
        stim_rec_dist = row.stim_rec_distance_center_to_center_mm;
        gain_value = row.gain_value;
        fs = 100000;
        dt = 1/fs;
        labchart_unit_conv = 1;
        stim_PW = row.stim_PW;
        stim_freq = 3.096 / 1000;
        V2mA = 0.1;
        Ntrig_expected = 10;
        phase_type = row.phase_type{1};
        t_start_sig = 0.6;
        
        % CV bounds (same as standard script)
        CV_bounds(1,1) = 5;   CV_bounds(1,2) = 100;
        CV_bounds_slowest = 0.1;
        
        y_min = -1300;
        y_max = -y_min;
        y_max_f = 10000;
        
        % Create output directory if it doesn't exist
        if ~exist(file_prefix_outputs, 'dir')
            mkdir(file_prefix_outputs);
        end
        
        % Check if cached raw data exists before attempting to process
        raw_cache_root = '/Volumes/bme_pena/Data_from_Duke_Collabs/CAP2/mat_cache_of_raw_data';
        first_channel_cache = fullfile(raw_cache_root, file_prefix_outputs(end-length(row.file_prefix_outputs_values{1})+1:end), sprintf('ch%d.mat', ch_num_neural_values(1)));
        
        if ~exist(first_channel_cache, 'file')
            fprintf('  SKIP: Cached raw data not found: %s\n', first_channel_cache);
            fprintf('        Run AnalysisScript_CAPs on Windows first to create the cache.\n');
            continue;
        end
        
        % Run artifact-removed preprocessing FOR EACH CHANNEL
        for ch_idx = 1:numel(ch_num_neural_values)
            ch_num_neural = ch_num_neural_values(ch_idx);
            
            % Check if output file already exists
            output_file = sprintf('%s/CAPData_ch%d.mat', file_prefix_outputs, ch_num_neural);
            
            if skip_if_exists && exist(output_file, 'file')
                fprintf('  SKIP: Channel %d already processed (file exists: %s)\n', ...
                    ch_num_neural, output_file);
                continue;
            end
            
            fprintf('  Running AnalysisScript_CAPs_ArtifactRemoved for channel %d...\n', ch_num_neural);
            
            AnalysisScript_CAPs_ArtifactRemoved(expdate, file_prefix_adicht_full, ...
                file_prefix_outputs, CAP_avg_type, subtractnoise, zerobaseline, ...
                digitalLPF, digitalLPF_fc_bandpass, digitaloffsetremoval, ch_num_Vinput, ...
                ch_num_neural, stim_rec_dist, gain_value, fs, dt, ...
                labchart_unit_conv, stim_PW, stim_freq, V2mA, Ntrig_expected, phase_type, ...
                t_start_sig, CV_bounds, CV_bounds_slowest, y_min, y_max, y_max_f, ...
                amp_values_nominal, amp_values_actual, trial_values);
        end
        
        fprintf('  Done with %s\n', file_prefix_outputs);
    else
        fprintf('  SKIP: TestPulseData file not found: %s\n', amp_file);
    end
end

fprintf('\nAll datasets processed.\n');
