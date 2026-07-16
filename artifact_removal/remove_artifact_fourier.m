function [u_time, x1_time, delta_t] = remove_artifact_fourier(c1, c2, fs)
    % Fourier-domain paired-polarity artifact removal
    % c1 = positive polarity recording
    % c2 = negative polarity recording
    % fs = sampling rate (Hz)
    
    %% --- Ensure column vectors ---
    c1 = c1(:);
    c2 = c2(:);
    N  = numel(c1);
    
    %% --- FFT of both recordings ---
    C1 = fft(c1);
    C2 = fft(c2);
    
    %% --- Estimate CAP time shift using cross-correlation ---
    [xcorr_val, lags] = xcorr(c2, c1);
    [~, idx] = max(abs(xcorr_val));
    delta_t = lags(idx) / fs;  % seconds?
    disp(delta_t)
    delta_t = .08;

    %% --- Visualize the CAP shift found by cross-correlation ---
    t = (0:N-1)/fs * 1000;   % ms
    
    % Shift c1 by delta_t (using interpolation)
    t_shifted = t - delta_t*1000;   % convert seconds → ms
    c1_shifted = interp1(t, c1, t_shifted, 'linear', 0);
    
    figure;
    subplot(2,1,1);
    plot(t, c1, 'k'); hold on;
    plot(t, c2, 'r');
    title(sprintf('Raw signals (delta_t = %.4f ms)', delta_t*1000));
    xlabel('Time (ms)'); ylabel('Amplitude');
    legend('c1 (positive)', 'c2 (negative)');
    
    subplot(2,1,2);
    plot(t, c2, 'r'); hold on;
    plot(t, c1_shifted, 'b');
    title('c1 shifted by delta_t to align with c2');
    xlabel('Time (ms)'); ylabel('Amplitude');
    legend('c2', 'c1 shifted');

    %% --- Frequency axis (MUST be row vector) ---
    omega = 2*pi*(0:N-1)/N;     % 1×N row vector
    
    %% --- Phase shift for each frequency ---
    theta = omega * delta_t;    % 1×N row vector
    
    %% --- Avoid division-by-zero when 1 + exp(jθ) = 0 ---
    den = 1 + exp(1j*theta);
    den(abs(den) < 1e-12) = 1e-12;   % regularization
    
    %% --- Artifact in frequency domain ---
    U = (exp(1j*theta) .* C1.' - C2.') ./ den;   % force row-vector math
    
    %% --- CAP (artifact-free) for polarity 1 ---
    X1 = (C1.' + C2.') ./ den;
    
    %% --- Back to time domain ---
    u_time  = real(ifft(U));
    x1_time = real(ifft(X1));
    
    %% --- Time axis ---
    t = (0:N-1)/fs * 1000; % ms
    
    %% --- Plot results ---
    figure;
    subplot(3,1,1);
    plot(t, c1, 'k'); hold on; plot(t, c2, 'r');
    title('Raw paired-polarity recordings');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
    subplot(3,1,2);
    plot(t, u_time, 'm');
    title('Estimated artifact template u(t)');
    xlabel('Time (ms)'); ylabel('Amplitude');
    
    subplot(3,1,3);
    plot(t, x1_time, 'b');
    title('Artifact-free CAP estimate x_1(t)');
    xlabel('Time (ms)'); ylabel('Amplitude');

end
