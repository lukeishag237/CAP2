%{

Return the difference between the largest positive peak and the largest
negative peak as the peak-to-peak value of the signal.

Alternatively (based on return_str value), return the average time point
between the largest positive peak and the largest negative peak.

%}
function return_value = find_CAP_amplitude_or_location(t,y,return_str)
Fs = 1/mode(diff(t));

% negative y to find the negative peak
min_prominence = 10;
[positive_peak_vals_all,positive_peak_loc_ms_all] = ...
    findpeaks(y,t,'MinPeakProminence',min_prominence);
[negative_peak_vals_all,negative_peak_loc_ms_all] = ...
    findpeaks(-y,t,'MinPeakProminence',min_prominence);

% calculate peak-to-peak value as the difference between the largest
% positive peak and the largest negative peak
[max_positive_peak_val,max_peak_ind] = max(positive_peak_vals_all);
[max_negative_peak_val,min_peak_ind] = min(-negative_peak_vals_all);
peak_to_peak = max_positive_peak_val-max_negative_peak_val;

% calculate the time of the AP as the mean of the largest
% positive peak and the largest negative peak
ap_latency_ms = mean([positive_peak_loc_ms_all(max_peak_ind),...
    negative_peak_loc_ms_all(min_peak_ind)]);

% plot gut check
flag_plot_gut_check = 0;
if (flag_plot_gut_check && strcmp(return_str,'peak_to_peak_value'))
    figure;
    plot(t,y);
    hold on;
    plot([ap_latency_ms,ap_latency_ms],[max_positive_peak_val,max_negative_peak_val]);
end

% return whatever specified by the return_str
if (strcmp(return_str,'peak_to_peak_value'))
    return_value = peak_to_peak;
elseif (strcmp(return_str,'inter_peak_to_peak_time'))
    return_value = ap_latency_ms;
elseif (strcmp(return_str,'negative_peak_time'))
    return_value = negative_peak_loc_ms_all(min_peak_ind);
else
    error('invalid value of return_str: %s',return_str)
end