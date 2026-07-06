%{

Use interpft to resample (upsample or downsample) the provided signal and time vector by the
specified upsample factor

original_signal must be either a column vector or a matrix
If original_signal is a matrix, this function operates on the columns

Below is a test example that should produce either an upsampled or
downsampled sinusoid depending on the chosen resampling factor (>1
upsamples; <1 downsamples)
    dt_original = 0.001;
    time_vector_original = 0:dt_original:1;
    signal_original = sin(2*pi*1*(time_vector_original));
    figure; plot(time_vector_original,signal_original,'.-')

    resample_factor = 0.05;
    [resampled_signal,resampled_time_vector] = ...
        resample_signal(signal_original',time_vector_original,...
        resample_factor);
    hold on; plot(resampled_time_vector,resampled_signal,'.')

%}
function [resampled_signal,resampled_time_vector] = resample_signal(original_signal,original_time_vector,resample_factor,type)

if (~exist('type','var') || isempty(type))
    type = 'ft';
end

if (strcmp(type,'spline'))
    dt_old = mode(diff(original_time_vector));
    dt_new = dt_old/resample_factor;
    resampled_time_vector = [original_time_vector(1):dt_new:original_time_vector(end)]';
    resampled_signal = interp1(original_time_vector,original_signal,resampled_time_vector,'spline');
else

assert(size(original_signal,1)~=1,'original_signal must be either a column vector or a matrix')

%%% Append a "taper" at the end of the signal so that the assumption of
%%% periodicity that interpft makes does not introduce high frequency artifacts
%%% in the resampled signal
% set the taper size such that you get very close to an integer value for
% the N resample (i.e., where N resample = (size(original_signal,1)+taper_size)*resample_factor )
taper_size_range = [0.25 0.5]; % if relative taper size is 0.5, then the taper will be roughly in the first harmonic; make it one of the low harmonics
N = size(original_signal,1); 
relative_taper_size_options = (1:N)/N; 
remainder_of_N_plus_taper_size_over_one = rem(((N+1):(2*N))*resample_factor,1); 
indices_of_valid_relative_taper_sizes = find(relative_taper_size_options>=taper_size_range(1) & relative_taper_size_options<=taper_size_range(end)); 
[~,min_idx] = min(remainder_of_N_plus_taper_size_over_one(indices_of_valid_relative_taper_sizes)) ; 
taper_size_relative_to_original_signal = relative_taper_size_options(indices_of_valid_relative_taper_sizes(min_idx)); %0.5; 
assert(taper_size_relative_to_original_signal<=taper_size_range(end));
taper_size = round(N*taper_size_relative_to_original_signal);
interp_X = [1:N,N+taper_size+(1:N)]; interp_Xq = (N+1):(N+taper_size);
taper = interp1(interp_X,[original_signal; original_signal],interp_Xq','spline'); %[1:size(original_signal,1), [1:n_repeats,n_repeats+taper_size+(1:n_repeats)],[zeros(1,n_repeats),ones(1,n_repeats)],n_repeats+(1:taper_size))'; % spline taper
extended_signal = [original_signal; taper]; % stitch the taper to the end of the signal
N_resampled = round(size(extended_signal,1)*resample_factor); %extended_signal_duration/(dt_original/resample_factor)

%%% Resample the signal
resampled_signal = interpft(extended_signal,N_resampled);

%%% Calculate time vector for resampled signal
ARBITRARY_PRECISION = 9; 
dt_original = round(mode(diff(original_time_vector)),ARBITRARY_PRECISION);
% extended_signal_time_vector = original_time_vector(1) + ((1:size(extended_signal,1))-1)*dt_original;
extended_signal_duration = size(extended_signal,1)*dt_original;
dt_resampled = extended_signal_duration/N_resampled;
resampled_time_vector = original_time_vector(1) + ((1:size(resampled_signal,1))-1)*dt_resampled;

%%% Discard the taper and any other extraneous time series data
resampled_signal(resampled_time_vector>original_time_vector(end) | resampled_time_vector<original_time_vector(1),:) = [];
resampled_time_vector(resampled_time_vector>original_time_vector(end) | resampled_time_vector<original_time_vector(1)) = [];
end

