function AnalysisScript_CAPs_Wrapper(expdate,file_prefix_outputs,CAP_avg_type,scaled_noise,flag_attenuation,...
                                                subtractnoise,zerobaseline,digitalLPF,digitalLPF_fc_bandpass,digitaloffsetremoval,ch_num_Vinput,...
                                                ch_num_neural_values,monopolarrec,experiment_metadata_table_row_i,raw_data_directory)

scrnsz = get(0,'screensize');

%% Define data targets
stim_rec_dist       = experiment_metadata_table_row_i.stim_rec_distance_center_to_center_mm;                   % [mm]
gain_value          = experiment_metadata_table_row_i.gain_value;
fs                  = 100000;               % [samples/s]
dt                  = 1/fs;
labchart_unit_conv  = 1;                    % to get [V]
stim_PW             = experiment_metadata_table_row_i.stim_PW;           % [ms]
stim_freq           = 3.096 / 1000;         % [kHz]
V2mA                = 0.1;                  % [mA/V]
Ntrig_expected      = 10;
phase_type          = experiment_metadata_table_row_i.phase_type{1};
t_start_sig         = 0.6;

CV_bounds(1,1) = 5;       CV_bounds(1,2) = 100; % [m/s]
% CV_bounds(2,1) = 1;         CV_bounds(2,2) = 1.5;
% CV_bounds(3,1) = 1.5;       CV_bounds(3,2) = 2;
% CV_bounds(4,1) = 0.5;       CV_bounds(4,2) = 2;
% CV_bounds(5,1) = 2;         CV_bounds(5,2) = 100;

CV_bounds_slowest   = 0.1;

trialnum_start_values   = experiment_metadata_table_row_i.trialnum_start_values;

file_prefix_metadata        = [raw_data_directory,'/',experiment_metadata_table_row_i.file_prefix_metadata{1}];
file_prefix_adicht          = [file_prefix_metadata, '/', experiment_metadata_table_row_i.file_prefix_adicht{1}];

y_min = -1300;
y_max = -y_min;
y_max_f = 10000;

% Attenuation
% amp_actual (mA) = p1 * amp_target + p2
% column 1: p1
% column 2: p2
attenuation_coeffs = [1 1];

%% Call analysis code
trialnum_start = trialnum_start_values;

% Load stimulation amplitudes & trial numbers, both already sorted in
% ascending order of stim amp
load([file_prefix_metadata 'TestPulseData_Trial' num2str(trialnum_start) '_ch' num2str(ch_num_neural_values(1)) '.mat'],'amp_values_nominal','trial_values')

if flag_attenuation == 0
    amp_values_actual = amp_values_nominal;
elseif flag_attenuation == 1
    p1                      = attenuation_coeffs(1);
    p2                      = attenuation_coeffs(2);
    amp_values_actual       = p1 .* amp_values_nominal + p2;
end

for ch_num_neural_ind = 1:numel(ch_num_neural_values)
    ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
    AnalysisScript_CAPs(expdate,file_prefix_adicht,file_prefix_outputs,...
                        CAP_avg_type,scaled_noise,subtractnoise,zerobaseline,...
                        digitalLPF,digitalLPF_fc_bandpass,digitaloffsetremoval,ch_num_Vinput,...
                        ch_num_neural,stim_rec_dist,gain_value,fs,dt,...
                        labchart_unit_conv,stim_PW,stim_freq,V2mA,Ntrig_expected,phase_type,...
                        t_start_sig,CV_bounds,CV_bounds_slowest,y_min,y_max,y_max_f,...
                        amp_values_nominal,amp_values_actual,trial_values)

    close all
end

%% If used bipolar or tripolar recording cuff but recorded from each contact 
% independently, then calculate CAPs, plot CAPs, plot recruitment curve,
% and fit recruitment curve

% ***** Load CAPs
for ch_num_neural_ind = 1:numel(ch_num_neural_values)
    ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
    fname = [file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'];
    load(fname,'ch_neural_avg')
    ch_neural_avg_all(ch_num_neural_ind,:,:) = ch_neural_avg;
end

% ***** Calculate bipolar or tripolar CAPs
if monopolarrec == 1
    if numel(ch_num_neural_values) == 2 % bipolar rec cuff
        ch_neural_avg_all(end+1,:,:) = ch_neural_avg_all(1,:,:) - ch_neural_avg_all(2,:,:);
        rec_type = 'bipolar';

    elseif numel(ch_num_neural_values) == 3 % tripolar rec cuff
        ch_neural_avg_all(end+1,:,:) = ch_neural_avg_all(2,:,:) - ((ch_neural_avg_all(1,:,:) + ch_neural_avg_all(3,:,:))/2);
        rec_type = 'tripolar';

    else
        error('need to write code for this case')
    end
end

% ***** Plot bipolar or tripolar CAPs
if monopolarrec == 1
    fname = [file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'];
    load(fname,'t_window','CAPvsTime_axlim','CV_latencybounds','y_CV_bounds','t_bounds_slowest','subplot_rows_all','subplot_cols_all')
    ch_num_neural = size(ch_neural_avg_all,1);

    figure
    set(gcf,'position',scrnsz)
    plot_CAPvsTime(numel(trial_values),squeeze(ch_neural_avg_all(end,:,:)),subplot_rows_all,subplot_cols_all,...
                            t_window,size(CV_bounds,1),CV_latencybounds,y_CV_bounds,...
                            t_bounds_slowest,CAPvsTime_axlim,expdate,ch_num_neural,...
                            phase_type,stim_PW,trial_values,amp_values_nominal)

    fname = [file_prefix_outputs '/CAPvsTime_allchannels'];
    saveas(gcf,[fname '.jpg']);
    saveas(gcf,[fname '.fig']);
end
                  
% ***** Overlay max CAPs
trial_ind = numel(trial_values);
figure; hold on

% Individual monopolar channels
for ch_num_neural_ind = 1:numel(ch_num_neural_values)
    ch_num_neural = ch_num_neural_values(ch_num_neural_ind);
    plot(t_window,ch_neural_avg_all(ch_num_neural_ind,:,trial_ind),'LineWidth',3,'DisplayName',['ch' num2str(ch_num_neural)])
end
axis(CAPvsTime_axlim)
axis square
title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
       [phase_type ', PW = ' num2str(stim_PW) ' \mus'],...
       ['Trial ' num2str(trial_values(trial_ind)) '; ' ...
       num2str(amp_values_nominal(trial_ind)) ' mA']})
xlabel('Time [ms]')
ylabel('Neural rec [\muV]')
legend

% Bipolar or tripolar max CAP
if monopolarrec == 1
    plot(t_window,ch_neural_avg_all(end,:,trial_ind),'LineWidth',3,'DisplayName',rec_type)
end

fname = [file_prefix_outputs '/CAPvsTime_maxCAP_allchannels'];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

% ***** Recruitment curve for post-hoc calculated bipolar or tripolar CAP
%{
NEED TO FILL THIS IN
%}

close all
