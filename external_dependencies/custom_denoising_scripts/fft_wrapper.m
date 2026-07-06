%{

Author: Edgar Pena
Date: Tue 12 Apr 2022

Wrapper to MATLAB's built-in fft function. It calculates a frequency vector
based on the provided fs. MATLAB's built-in fft function is agnostic to
this information, so it doesn't calculate freq, yet freq is
very useful for interpreting fft results.

TODO before putting this on MATLAB Central:
- check if others did something similar
- maybe expand to calculate Fourier Series coefficients too

%}

function [Y,freq] = fft_wrapper(X, N, DIM, fs)

if ~exist('DIM','var') || isempty(DIM)
    DIM = 1;
end

if ~exist('N','var') || isempty(N)
    N = size(X,DIM);
end

%%% Run FFT using MATLAB's built-in fft function and built-in inputs
Y = fft(X,N,DIM);

%%% Calculate freq using fs
df = fs/N;
if (mod(N,2)==0)
    % for even FFT, the first half starts with the 0 Hz component, while the
    % second half starts with a real (i.e., not complex-valued)  high frequency
    % component; every other value has a complex conjugate pair
    freq = [((1:(N/2))-1)*df, ((1:(N/2)))*df - (N/2 + 1)*df];
    
else
    % for odd FFT, the first half starts with the 0 Hz component, and the last
    % value of the first half is the complex conjugate of the first value of
    % the second half
    freq = [((1:ceil(N/2))-1)*df, ((1:floor(N/2)))*df - ceil(N/2)*df];
end
