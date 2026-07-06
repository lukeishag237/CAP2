%{

Description: 
'windowed_signal = split_signal(signal,num_samples_per_window,num_samples_overlap))' returns 
 a windowed version of the input signal, 'windowed_signal'. The signal is
 split into segments of length 'num_samples_per_window', with the amount of
 overlap specified by 'num_samples_overlap'. Each column of 'windowed_signal' 
is a separate snippet of the signal.

Input Arguments:
signal (vector) - Input signal, specified as either a row vector or column vector
num_samples_per_window (scalar) - indicates the number of samples per window
num_samples_overlap (scalar) (optional)  - Indicates the number of samples of overlap

Output Arguments:
windowed_signal - windowed signal. Each column is a separate snippet of the signal
index_matrix - matrix defined such that windows_signal = signal(index_matrix)

%}

function [windowed_signal,index_matrix] = split_signal(signal,num_samples_per_window,num_samples_overlap)

assert(nargin>=2,'Error: Must at least specify signal and size of windows')

if (nargin<3)
    num_samples_overlap = 0;
end

%%% Specify start and end indices
start_indices = 1:(num_samples_per_window-num_samples_overlap):length(signal);
end_indices = start_indices+num_samples_per_window-1;

%%% Remove invalid indices
invalid_indices = find(end_indices>length(signal));
start_indices(invalid_indices) = [];
end_indices(invalid_indices) = [];

%%% Expand all indices (M-by-N matrix, where M is num_samples_per_window,
%%% and N is the number of windows produced from the specified signal and
%%% overlaps
index_matrix = cell2mat(arrayfun(@(x,y) (x:y)', start_indices, end_indices,'UniformOutput',false));

%%% Create windowed signal
windowed_signal = signal(index_matrix);

