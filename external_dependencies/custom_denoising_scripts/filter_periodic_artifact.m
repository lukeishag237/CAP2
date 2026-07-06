%{

Filters periodic artifacts by estimating the artifact on a period-by-period
basis. It does this by resampling the signal using the true artifact
period, then estimating the artifact at each window using adjacent windows,
and finally subtracting the artifact out

data_tmp - signal to be filtered
fs - sampling rate (in Hz)
nominal_noise_frequency - nominal frequency of the noise

%}

function [data_tmp_filtered,time_vector_sec] = filter_periodic_artifact(data_tmp,fs_Hz,nominal_noise_frequency)
%%

dt = 1/fs_Hz; % [sec]

time_vector_sec = ((1:length(data_tmp))-1)*dt; % [sec]

%%% Start preprocessing signal for artifact removal
input_signal = data_tmp;

flag_determine_precise_60Hz_noise_freq = 1;
if (flag_determine_precise_60Hz_noise_freq)
    %%% Figure out the main 60 Hz frequency precisely
    fft_factor = 80;
    nfft = length(input_signal)*fft_factor;
    [fft_data_tmp,freq_vector] = fft_wrapper(input_signal,nfft,[],fs_Hz);
    EXPECTED_NOISE_RANGE = nominal_noise_frequency + [-1 1]; % [Hz]
    main_60Hz_noise_range_indices = find(freq_vector>EXPECTED_NOISE_RANGE(1) & freq_vector<EXPECTED_NOISE_RANGE(2));
    [~,max_idx] = max(abs(fft_data_tmp((main_60Hz_noise_range_indices))));
    main_60_noise_freq_Hz = freq_vector(main_60Hz_noise_range_indices(max_idx)); % [Hz]
    clearvars fft_data_tmp freq_vector
else
    main_60_noise_freq_Hz = 60; % [Hz]
end

flag_resample_signal = 1;
if (flag_resample_signal)
    %%% Resample the signal such that the period of the main 60 Hz
    %%% noise is an integer multiple of the resampled dt
    resample_factor = round(((1/main_60_noise_freq_Hz)/(1/fs_Hz)))/((1/main_60_noise_freq_Hz)/(1/fs_Hz));
    [resampled_signal,resampled_time_vector] = resample_signal(input_signal,time_vector_sec,resample_factor);
    resampled_dt = mode(diff(resampled_time_vector));
else
    resample_factor = 1;
    resampled_signal = input_signal;
    resampled_time_vector = time_vector_sec;
    resampled_dt = mode(diff(resampled_time_vector));
end

%%% Run LPF to remove 60 Hz noise
% As of Tue 12 Apr 2022, filter_method 4 seems the best among these
% options. Try running the lines below after the post-processing step
% to verify visually. The simplest possible option that denoises is the
% FIR, although it also smooths out sharp transitions, and isn't
% particularly fast. (An IIR alternative would likely be faster than
% the FIR method of filter_method 3.)
%{
    toc; figure('units','normalized','outerposition',[0 0 1 1]);
    plot(data_tmp); hold on; plot(data_tmp_filtered);  title(sprintf('case = %d',filter_method));
    
%}
filter_method = 4;
switch filter_method
    case 0 % none
        resampled_data_tmp_filtered = resampled_signal;
        resampled_data_tmp_filtered_time_vector = ((1:length(resampled_data_tmp_filtered))-1)*resampled_dt;
        
    case 1 % FIR filter
        b = fir1(10e3,60/(resample_factor*fs_Hz/2),'low');
        resampled_data_tmp_filtered = filtfilt(b,1,resampled_signal);
        resampled_data_tmp_filtered_time_vector = ((1:length(resampled_data_tmp_filtered))-1)*resampled_dt;
        
    case 2 % just 'jiggle' the signal and average the 'jiggled' signals
        n = 64;
        resampled_data_tmp_filtered = zeros(size(resampled_signal));
        for n_val_i = -n:1:n
            resampled_data_tmp_filtered = resampled_data_tmp_filtered + resampled_signal(circshift(1:length(resampled_signal),n_val_i));
        end
        resampled_data_tmp_filtered = resampled_data_tmp_filtered/(2*n+1);
        resampled_data_tmp_filtered_time_vector = ((1:length(resampled_data_tmp_filtered))-1)*resampled_dt;
        
    case 3 % an adaptive template subtraction: fit the noise signal with a fourier fit (adjusting its timing based on the data) then subtract the noise template from the data
        %%% Split the signal into windows of period equal to a single
        %%% cycle of 60 Hz noise, then remove the mean
        
        num_samples_per_window = round((1/main_60_noise_freq_Hz)/resampled_dt);
        num_samples_overlap = 0;
        [windowed_signal,~] = split_signal(resampled_signal,num_samples_per_window,num_samples_overlap);
        demeaned_windowed_signal = bsxfun(@minus,windowed_signal,mean(windowed_signal));
        
        %%% Use the first window to approximate the artifact shape
        FIRST_WINDOW_IDX = 1;
        [fft_window_i,freq_vector] = fft_wrapper(demeaned_windowed_signal(:,FIRST_WINDOW_IDX),[],[],1/resampled_dt);
        harmonic_limit_Hz = 750; % [Hz] Include up to 720 Hz harmonics of the 60 Hz artifact (so keep up to 750 Hz)
        fft_window_i_nulled = fft_window_i; fft_window_i_nulled(abs(freq_vector)>harmonic_limit_Hz)=0;
        artifact_template = ifft(fft_window_i_nulled)';
        
        %%% Make matrix with all the shited templates (only shift up to +- 20)
        shift_amount = 20;
        x = circshift(artifact_template,-shift_amount);
        r = x;
        c = [x(1) fliplr(x(end-size(demeaned_windowed_signal,1)+2:end))];
        xConv = toeplitz(c(1:(2*shift_amount)),r);
        
        %%% Multiply shifted templates to de-meaned windowed signals
        gg=xConv*demeaned_windowed_signal;
        [~,max_idx] = max(gg);
        smooth_max_idx = round(movmedian(max_idx,15)); % smooth the shifts
        indices_matrix = -(smooth_max_idx-(shift_amount+1))'+(1:length(artifact_template));
        indices_matrix(indices_matrix>length(artifact_template)) = indices_matrix(indices_matrix>length(artifact_template)) - length(artifact_template);
        indices_matrix(indices_matrix<=0) = indices_matrix(indices_matrix<=0) + length(artifact_template);
        artifact_template_matrix = artifact_template(indices_matrix)';
        
        %%% Subtract noise
        cleaned_windowed_signal = windowed_signal - artifact_template_matrix;
        
        %%% Reconstruct full signal
        resampled_data_tmp_filtered = reshape(cleaned_windowed_signal,[],1);
        resampled_data_tmp_filtered_time_vector = (1:length(resampled_data_tmp_filtered)-1)*resampled_dt;
        
    case 4 % another adaptive template subtraction: estimate the noise signals then subtract them from the original signals
        %%% From the resampled signal, calculate the FFT and use only
        %%% the high frequencies (i.e., greater than a specified cutoff
        %%% frequency) to generate a noise signal in the time domain.
        %%% Then use the noise signal and a moving median window to
        %%% calculate artifacts for each snippet of the signal
        
        % Create a noise signal that only contains the high frequencies
        % (above cutoff) from the resampled signal
        resampled_dt = mode(diff(resampled_time_vector)); % [sec]
        [fft_resampled_signal,freq_vector] = fft_wrapper(resampled_signal,[],[],1/resampled_dt);
        offset_Hz = 6; % [Hz] amount below nominal noise frequency below which the target noise is not present at all
        cutoff_freq = nominal_noise_frequency - offset_Hz; % [Hz] zero out all frequency components greater than cutoff_freq
        fft_resampled_signal(abs(freq_vector)<=cutoff_freq) = 0;
        resampled_noise_signal = ifft(fft_resampled_signal);
        
        % Smooth the noise signal by taking the median of adjacent 60
        % Hz noise segments
        num_samples_per_window = round((1/main_60_noise_freq_Hz)/resampled_dt);
        num_samples_overlap = 0;
        [windowed_noise_signal,~] = split_signal(resampled_noise_signal,num_samples_per_window,num_samples_overlap);
        median_n_templates = 30; % each template will be the median of the present signal and the adjacent median_n_templates
        templates = movmedian(windowed_noise_signal',median_n_templates)';
        smoothed_resampled_nose_signal = reshape(templates,[],1);
        
        % Subtract the smoothed resampled noise signal from resampled
        % signal, then resample back to the original
        resampled_data_tmp_filtered = resampled_signal(1:length(smoothed_resampled_nose_signal))-smoothed_resampled_nose_signal;
        resampled_data_tmp_filtered_time_vector = ((1:length(resampled_data_tmp_filtered))-1)*resampled_dt;
        
        
    otherwise
        error('invalid switch case value')
end

%%% Post-process the signal for artifact removal
% Resample to original temporal resolution
if (flag_resample_signal)
    [data_tmp_filtered,time_vector_sec] = resample_signal(resampled_data_tmp_filtered,resampled_data_tmp_filtered_time_vector,1/resample_factor);
else
    data_tmp_filtered = resampled_data_tmp_filtered;
    time_vector_sec = resampled_data_tmp_filtered_time_vector;
end

flag_moving_average_smoothing = 1;
if (flag_moving_average_smoothing)
    % Smooth the signal out with a moving average (choose moving
    % average instead of moving median because the moving average
    % smoothed out the random jumps in the signal a bit better; but
    % both are very similar) A moving window of size 10 achieved
    % visually better smoothing than a window size of 2 or 4, while
    % still maintaining the main features of the signal (based on
    % visual inspection)
    moving_window_size = 10; % number of samples
    data_tmp_filtered = movmean(data_tmp_filtered,moving_window_size);
end
