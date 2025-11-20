function AnalysisScript_CAPs(expdate,file_prefix_adicht,file_prefix_outputs,...
                    CAP_avg_type,scaled_noise,subtractnoise,zerobaseline,...
                    digitalLPF,digitalLPF_fc_bandpass,digitaloffsetremoval,ch_num_Vinput,...
                    ch_num_neural,stim_rec_dist,gain_value,fs,dt,...
                    labchart_unit_conv,stim_PW,stim_freq,V2mA,Ntrig_expected,phase_type,...
                    t_start_sig,CV_bounds,CV_bounds_slowest,y_min,y_max,y_max_f,...
                    amp_values_nominal,amp_values_actual,trial_values)

scrnsz      = get(0,'screensize');
indiv_plots = 0;

%% Define data targets
% tmp                 = sort(trial_values);
% trial_num_start     = tmp(1);

Ntrials             = numel(trial_values);
N_CVranges          = size(CV_bounds,1);

% Voltage input to voltage-to-current converter [V]
Vinput_values       = amp_values_actual./V2mA;    

subplot_rows        = 3;
subplot_cols        = 1;

% Subplot dimensions to plot all trials for given block in one figure
subplot_rows_all = 3;
subplot_cols_all = 7;

check_Nsubplots = 1;
while (check_Nsubplots == 1)
    if ((subplot_rows_all*subplot_cols_all) < Ntrials)
        subplot_rows_all = subplot_rows_all + 1;
    else
        check_Nsubplots = 0;
    end
end

%% Time and CV parameters
t_bounds_slowest    = stim_rec_dist./CV_bounds_slowest;

% Latency bounds
CV_latencybounds        = stim_rec_dist./CV_bounds;


% % Check that fastest CV bound isn't evidently in stim artifact
% if strcmp(phase_type,'monophasic') == 1
%     if min(CV_latencybounds(:)) < stim_PW
%         error('Looking for signal within stim artifact')
%     end
% elseif strcmp(phase_type,'biphasic') == 1
%     if min(CV_latencybounds(:)) < 2*stim_PW
%         error('Looking for signal within stim artifact')
%     end
% else
%     error('need to write code for this case')
% end

% The start of the noise window should be later than the expected end of
% the C-type CAP
stim_period             = 1/stim_freq;          % [ms]
t_max                   = t_bounds_slowest;     % [ms]
earliest_tstart_noise   = stim_period - t_max;

pretime       = earliest_tstart_noise;    % Time before event to include in window [ms]
                                                    % Negative to specify window starting after event
posttime      = t_max;                         % Time after event to include in window [ms]

% t_window
pretime_idx     = round(pretime*fs/1000);
posttime_idx    = round(posttime*fs/1000);
Nidx            = numel(-pretime_idx:posttime_idx);
t_window        = -dt*pretime_idx:dt:dt*posttime_idx; t_window = t_window.*1e3;
CV_window       = stim_rec_dist ./ t_window;

% Convert latency bounds to indices
CV_latencybounds_idx    = zeros(size(CV_latencybounds));
for CVrange_ind = 1:N_CVranges
    for bound_ind = 1:2 % Lower and upper bounds on each range
        [~,CV_latencybounds_idx(CVrange_ind,bound_ind)] = min(abs(t_window-CV_latencybounds(CVrange_ind,bound_ind)));
    end
end

% Latency duration for each range of CVs
CV_latencydurations = -diff(CV_latencybounds,1,2);

% Noise for each trial
clear noise_latencybounds
if scaled_noise == 1
    noise_latencybounds(1) = -earliest_tstart_noise; %-39; % -39 ms was evaluated for suitability in re-analysis of Pelot & Grill 2020
    noise_latencybounds(2) = -1;  % Avoid stim artifact in case rising edge (at t=0) is second phase
    
elseif scaled_noise == 0
    noise_latencybounds(1) = -earliest_tstart_noise; %-150;
    noise_latencybounds(2) = -1;  % NOT USED
end

noise_durn           = noise_latencybounds(2) - noise_latencybounds(1);
noise_multipliers    = CV_latencydurations ./ noise_durn;

if subtractnoise == 1
    % ***** Check validity of noise bounds
    
    % Noise bounds should be negative, and more specifically, less than -PW
    % (i.e., before expected stim artifact (if rising edge has negative pulse
    % before t = 0)
    if (noise_latencybounds(1) > -stim_PW) || (noise_latencybounds(2) > -stim_PW)
        error('Invalid noise time bounds; must be <-PW to avoid overlap with stim artifact')
    end

    % Noise bounds should be greater than pretime, i.e., within the calculated CAP
    if (noise_latencybounds(1) < -pretime) || (noise_latencybounds(2) < -pretime)
        error('Invalid noise time bounds; must be within calculated CAP')
    end

    % If NOT approximating (i.e., scaling) the noise, then the noise bounds
    % must leave enough time for the longest latency range of the signals
    if scaled_noise == 0
        if noise_durn > max(CV_latencydurations)
            error('Invalid noise time bounds; available noise_durn is shorter than the longest latency range defined by CV_bounds; switch to scaled_noise = 1')
        end
    end

    % The start of the noise window should be later than the expected end of
    % the C-type CAP
    if (noise_latencybounds(1) < -earliest_tstart_noise)
        error('Invalid noise time bounds; must avoid overlapping noise with CAP')
    end
end

% Convert noise bounds from milliseconds to indices
noise_latencybounds_idx = zeros(size(noise_latencybounds));
for bound_ind = 1:2 % Lower and upper bounds on each range
    [~,noise_latencybounds_idx(bound_ind)] = min(abs(t_window-noise_latencybounds(bound_ind)));
end

%% Plotting parameters
thresh_targets  = 0.5;%[0.2 0.5 0.8];
cmp_thresh      = parula;
Ncmp            = numel(thresh_targets);
cmp_thresh      = cmp_thresh(1:floor(size(cmp_thresh,1)/Ncmp):end,:);

% Subplot dimensions for a panel for each CV range
subplot_rows_CVranges = 2;
subplot_cols_CVranges = 3;
if ((subplot_rows_CVranges*subplot_cols_CVranges) < N_CVranges)
    error('Need more subplots for N_CVranges')
end

cmp_CVranges    = parula;
Ncmp            = N_CVranges;
cmp_CVranges    = cmp_CVranges(1:floor(size(cmp_CVranges,1)/Ncmp):end,:);

y_CV_bounds = y_min:1:y_max;

% Axis limits
CAPvsTime_axlim = [-0.5 max(CV_latencybounds(:)) y_min y_max];
CAPvsCV_axlim   = [0    max(CV_bounds(:))        y_min y_max];
CAPvsFreq_axlim = [0    4000                     0     y_max_f];

%% Load data
fname       = file_prefix_adicht;
mydata      = adi.readFile(fname);
records     = mydata.records;       % New record every time LabChart is started/stopped
comments    = [records.comments];   % Concatenate all comments from all records

ch_Vinput_all = mydata.channel_specs(ch_num_Vinput);
ch_neural_all = mydata.channel_specs(ch_num_neural);

% Contents and relevant metadata of comments
comment_text    = {comments.str};
comment_record  = [comments.record];
comment_time    = [comments.time];
Ncomments       = numel(comment_time);

%% Parse data to pull out blocks of data for each trial
trial_start_prefix = 'Trial ';
trial_end_prefix   = 'Trial ';
trial_end_suffix   = ' End';

trial_records = zeros([numel(trial_values) 1]);
trial_time = zeros([numel(trial_values) 2]);

for trial_ind = 1:numel(trial_values)
    comment_trial_start = [trial_start_prefix num2str(trial_values(trial_ind))];
    comment_trial_end   = [trial_end_prefix   num2str(trial_values(trial_ind)) trial_end_suffix];
    
    for comment_ind = 1:Ncomments
        if strcmp(comment_text{comment_ind},comment_trial_start)
            trial_records(trial_ind) = comment_record(comment_ind);
            trial_time(trial_ind,1)  = comment_time(comment_ind);
        elseif strcmp(comment_text{comment_ind},comment_trial_end)
            trial_time(trial_ind,2)  = comment_time(comment_ind);
        end
    end
    
    ch_Vinput{trial_ind} = ch_Vinput_all.getData(trial_records(trial_ind),'time_range',[trial_time(trial_ind,1) trial_time(trial_ind,2)]);
    ch_neural{trial_ind} = ch_neural_all.getData(trial_records(trial_ind),'time_range',[trial_time(trial_ind,1) trial_time(trial_ind,2)]);
    
    % Convert to V
%     ch_neural{trial_ind} = ch_neural{trial_ind} .* labchart_unit_conv;
    
    if indiv_plots == 1
        figure(trial_ind)
        set(gcf,'position',[scrnsz(3)*0.25 scrnsz(4)*0.1 scrnsz(3)*0.5 scrnsz(4)*0.8])
        subplot(subplot_rows,subplot_cols,1); plot(ch_Vinput{trial_ind}); title(comment_trial_start)
        subplot(subplot_rows,subplot_cols,2); plot(ch_neural{trial_ind});
    end
end

%% Stim-triggered averaging of neural data
threshold_def = 0.001; % Low noise in signal from NIDAQ, so low threshold [V]

ch_neural_avg   = zeros(Nidx,numel(trial_values));
Ntrig           = zeros(numel(trial_values),1);

for trial_ind = 1:numel(trial_values)
    % Collect indices of rising edges of Vinput (also falling edges if it
    % finds <Ntrig_expected)
    threshold = max([threshold_def Vinput_values(trial_ind)/10]);
%     [trig_indices] = thresh_index_SDcurve(ch_Vinput{trial_ind},threshold,Ntrig_expected);
    [~,trig_indices] = findpeaks(diff(abs(ch_Vinput{trial_ind})),'MinPeakHeight',0.01,'MinPeakDistance',0.9*(stim_period/1e3)*fs);
    
    if (strcmp(phase_type,'biphasic_burst') == 1)
        trig_indices = trig_indices(1:N_pulses_in_burst:end);
    end
    
    trig_indices_start = trig_indices - pretime_idx;
    trig_indices_end   = trig_indices + posttime_idx;

    % Exclude first trig_indices if it's too close to the beginning of the epoch
    if trig_indices_start(1) <= 0
        trig_indices_start(1) = [];
        trig_indices_end(1) = [];
    end
    
    % Exclude last trig_indices if it's too close to the end of the epoch
    if trig_indices_end(end) > numel(ch_Vinput{trial_ind})
        trig_indices_start(end) = [];
        trig_indices_end(end) = [];
    end

    % Raw data + LPF (LPF to mimic filtering KHF artifact)
    data_tmp = ch_neural{trial_ind};
    data_tmp = data_tmp .* labchart_unit_conv;  % Convert to V
    data_tmp = data_tmp .* 1e6;                 % Convert to uV
    if digitalLPF == 1
        if fs == 20000
            data_tmp = FilterDesigner_IIR_Chebyshev1_20ksamples_2kto5k(data_tmp);
        elseif fs == 100000
            data_tmp = FilterDesigner_IIR_Chebyshev1_100ksamples_2kto5k(data_tmp);
        elseif fs == 200000
            if digitalLPF_fc_bandpass == 2
                data_tmp = FilterDesigner_IIR_Chebyshev1_200ksamples_2kto5k(data_tmp);
            elseif digitalLPF_fc_bandpass == 15
                data_tmp = FilterDesigner_IIR_Chebyshev1_200ksamples_15kto18k(data_tmp);
            else
                error('Digital LPF not available for fs = 200 kHz with current digitalLPF_fc_bandpass')
            end
        else
            error('Digital LPF not available with current sampling freq')
        end
%     else
%         data_tmp_lpf = data_tmp;
    end
    
    % Average window of neural data around each value in trig_indices
    Ntrig(trial_ind) = numel(trig_indices_start);
    if Ntrig(trial_ind) ~= Ntrig_expected
        error('Wrong number of pulses detected')
    end

    clear tmp_pre
    for trig_ind = 1:Ntrig(trial_ind)
        start_ind = trig_indices_start(trig_ind);
        stop_ind  = trig_indices_end(trig_ind);
%         ch_neural_avg(:,trial_ind) = ch_neural_avg(:,trial_ind) + data_tmp(start_ind:stop_ind,1);
        tmp_pre(trig_ind,:) = data_tmp(start_ind:stop_ind,1);
    end
    
    if CAP_avg_type == 1        % Mean
        ch_neural_avg(:,trial_ind) = mean(tmp_pre,1);
%         ch_neural_avg(:,trial_ind) = ch_neural_avg(:,trial_ind)./Ntrig(trial_ind);
    elseif CAP_avg_type == 2    % Median
        ch_neural_avg(:,trial_ind) = median(tmp_pre,1);
    end

    ch_neural_avg(:,trial_ind) = ch_neural_avg(:,trial_ind)./gain_value;

    if indiv_plots == 1
        figure(trial_ind); subplot(subplot_rows,subplot_cols,3)
        plot(t_window,ch_neural_avg(:,trial_ind),'k')
        xlabel('Time [ms]')
        ylabel('Neural rec [\muV]')
        xlim([t_window(1) t_window(end)])
        %     tmp = ylim; tmp = max(abs(tmp)); ylim([-tmp tmp])
        ylim([y_min y_max])
    end
end

% disp('Ntrig: '); disp(num2str(Ntrig));
% if ((min(Ntrig) < Ntrig_expected-3) || (max(Ntrig) > Ntrig_expected))


clear ch_neural ch_Vinput

%% Remove baseline offset

for trial_ind = 1:Ntrials
    y_tmp       = ch_neural_avg(:,trial_ind);
    if zerobaseline == 0
        baseline_tmp = 0;
        
    elseif zerobaseline == 1
        % Why? Because recording force data with AC-coupling and no
        % HPF, so want to bring everything up to zero based on the
        % "quiet" time in the pre-KHF twitches
        ind_vector_noise_tmp   = noise_latencybounds_idx(1):noise_latencybounds_idx(2);
        baseline_tmp           = mean(y_tmp(ind_vector_noise_tmp));
    end
    
    ch_neural_avg(:,trial_ind) = y_tmp - baseline_tmp;
end

%% Plot stim-triggered averages vs. t
figure
set(gcf,'position',scrnsz)

plot_CAPvsTime(Ntrials,ch_neural_avg,subplot_rows_all,subplot_cols_all,...
                t_window,N_CVranges,CV_latencybounds,y_CV_bounds,...
                t_bounds_slowest,CAPvsTime_axlim,expdate,ch_num_neural,...
                phase_type,stim_PW,trial_values,amp_values_nominal)

fname = [file_prefix_outputs '/CAPvsTime_ch' num2str(ch_num_neural)];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

%% Plot stim-triggered averages vs. CV
figure
set(gcf,'position',scrnsz)

plot_CAPvsCV(Ntrials,ch_neural_avg,subplot_rows_all,subplot_cols_all,...
                CV_window,N_CVranges,CV_bounds,y_CV_bounds,...
                CV_bounds_slowest,CAPvsCV_axlim,expdate,ch_num_neural,...
                phase_type,stim_PW,trial_values,amp_values_nominal)

fname = [file_prefix_outputs '/CAPvsCV_ch' num2str(ch_num_neural)];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

%% Plot FFT of stim-triggered averages
%{
% ***** SEE NOTES IN LAB MEETING, 20220617
h1 = figure;
set(gcf,'position',scrnsz)

h2 = figure;
set(gcf,'position',scrnsz)

t_window_postartifact   = t_window(t_window > t_start_sig);
Npts                    = numel(t_window_postartifact);
f_window                = fs*(0:(Npts/2))/Npts;
    
for trial_ind = 1:Ntrials
    % Pull CAP data
    y_tmp_t = ch_neural_avg(:,trial_ind);
    y_tmp_t = y_tmp_t(end-numel(t_window_postartifact)+1:end);
    y_tmp_f = abs(fft(y_tmp_t));
    y_tmp_f = y_tmp_f(1:ceil(numel(y_tmp_f)/2));
    
    ch_neural_avg_f(trial_ind,:) = y_tmp_f;

    % ***** CAP vs. time without stim artifact
    set(0, 'CurrentFigure', h1)
    
    % Identify the correct subplot for the current trial
    subplot(subplot_rows_all,subplot_cols_all,trial_ind); hold on
    
    % Plot data
    plot(t_window_postartifact,y_tmp_t,'k')
    axis(CAPvsTime_axlim)
    
    if (trial_ind == 1)
        title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
               [phase_type ', PW = ' num2str(stim_PW) ' \mus'],...
               ['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA']})
        xlabel('Time [ms]')
        ylabel('Neural rec [\muV]')
    else
        title(['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA'])
    end
    
    % ***** CAP vs. freq
    set(0, 'CurrentFigure', h2)
    
    % Identify the correct subplot for the current trial
    subplot(subplot_rows_all,subplot_cols_all,trial_ind); hold on
    
    % Plot data
    plot(f_window,y_tmp_f,'k')
    axis(CAPvsFreq_axlim)
    
    if (trial_ind == 1)
        title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
               [phase_type ', PW = ' num2str(stim_PW) ' \mus'],...
               ['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA']})
        xlabel('Freq [Hz]')
        ylabel('FFT of Neural rec [\muV]')
    else
        title(['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA'])
    end
end 

fname = [file_prefix_outputs '/CAPvsTime_AfterArtifact_ch' num2str(ch_num_neural)];
saveas(h2,[fname '.jpg']);
saveas(h2,[fname '.fig']);

fname = [file_prefix_outputs '/CAPvsFreq_ch' num2str(ch_num_neural)];
saveas(h2,[fname '.jpg']);
saveas(h2,[fname '.fig']);
%} 

%% Calculate & plot AUC recruitment curves

clear AUC AUC_noise

figure
set(gcf,'position',[scrnsz(3)*0.25 scrnsz(4)*0.1 scrnsz(3)*0.6 scrnsz(4)*0.8])
        
for trial_ind = 1:Ntrials
    y_tmp       = ch_neural_avg(:,trial_ind);
    
    for CVrange_ind = 1:N_CVranges
        ind_vector_signal_tmp                       = CV_latencybounds_idx(CVrange_ind,2):CV_latencybounds_idx(CVrange_ind,1);
        AUC(CVrange_ind,trial_ind)                  = trapz(abs(y_tmp(ind_vector_signal_tmp)))*dt*1000;
        
        if subtractnoise == 0
            AUC_noise(CVrange_ind,trial_ind) = 0;

        elseif subtractnoise == 1
            if scaled_noise == 1
                ind_vector_noise_tmp                = noise_latencybounds_idx(1):noise_latencybounds_idx(2);
                noise_base                          = trapz(abs(y_tmp(ind_vector_noise_tmp)))*dt*1000;
                AUC_noise(CVrange_ind,trial_ind)    = noise_base * noise_multipliers(CVrange_ind);

            elseif scaled_noise == 0
                ind_vector_noise_tmp                = noise_latencybounds_idx(1,1):(noise_latencybounds_idx(1,1)+numel(ind_vector_signal_tmp)-1);
                AUC_noise(CVrange_ind,trial_ind)    = trapz(abs(y_tmp(ind_vector_noise_tmp)))*dt*1000;
            end
        end
    end

    AUC(:,trial_ind) = AUC(:,trial_ind) - AUC_noise(:,trial_ind);
    
    % One panel for each CV range
    for CVrange_ind = 1:N_CVranges
        subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind); hold on
        plot(amp_values_actual(trial_ind),AUC(CVrange_ind,trial_ind),'*','Color',cmp_CVranges(CVrange_ind,:))
        
        if (trial_ind == 1)
            title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
                [phase_type ', PW = ' num2str(stim_PW) ' \mus'],...
                [num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s']})
            xlabel('Stim amp (mA)')
            ylabel('AUC (\muV*s)')
        end
    end
    
    % Overlay recruitment curves for all CV ranges
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,subplot_rows_CVranges*subplot_cols_CVranges); hold on
    for CVrange_ind = 1:N_CVranges
        if (trial_ind == 1)
            plot(amp_values_actual(trial_ind),AUC(CVrange_ind,trial_ind),'*','Color',cmp_CVranges(CVrange_ind,:),...
                 'DisplayName',[num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s'])
        else
            plot(amp_values_actual(trial_ind),AUC(CVrange_ind,trial_ind),'*','Color',cmp_CVranges(CVrange_ind,:),...
                 'HandleVisibility','off')
        end
    end
    title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
           [phase_type ', PW = ' num2str(stim_PW) ' \mus']})
    xlabel('Stim amp (mA)')
    ylabel('AUC (\muV*s)')
    legend('Location','northwest')
end

%% Fit recruitment curves
rsquare_precision   = 3;
ft                  = fittype('fitfunc_recruitmentcurve(x,a,b,c)','independent','x','dependent','y');
opts                = fitoptions('Method','NonlinearLeastSquares');
opts.Display        = 'Off';

max_stim_amp        = max(amp_values_actual);
x_values_stimamp    = 0:max_stim_amp/100:max_stim_amp;

ft_AUC = []; % So that there's still a variable to save even if there aren't enough points for the fit

x_tmp   = amp_values_actual;

clear y_tmp y_tmp_max
for CVrange_ind = 1:N_CVranges
    y_tmp(CVrange_ind,:)    = AUC(CVrange_ind,:);
end

Nmin = 3; % Min number of data points to fit the equation
if numel(x_tmp) < Nmin
    disp('Too few points for sigmoid fit')
    
else
    for CVrange_ind = 1:N_CVranges
        opts.StartPoint     = [max(y_tmp(CVrange_ind,:))  20              max(x_tmp)/2];
        opts.Lower          = [max(y_tmp(CVrange_ind,:))/2  1               0];
        % opts.Upper          = [1000 200 Inf];
        
        subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)
        [ft_AUC{CVrange_ind}, gof_AUC{CVrange_ind}] = fit(x_tmp',y_tmp(CVrange_ind,:)',ft,opts);
        plot(x_values_stimamp,fitfunc_recruitmentcurve(x_values_stimamp,ft_AUC{CVrange_ind}.a,...
            ft_AUC{CVrange_ind}.b,ft_AUC{CVrange_ind}.c),'k--')
        disp([num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s'])
        ft_AUC{CVrange_ind}
        gof_AUC{CVrange_ind}
        xlim([0 max(x_tmp)*1.1])
    end
    
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,subplot_rows_CVranges*subplot_cols_CVranges); hold on
    for CVrange_ind = 1:N_CVranges
        plot(x_values_stimamp,fitfunc_recruitmentcurve(x_values_stimamp,ft_AUC{CVrange_ind}.a,...
            ft_AUC{CVrange_ind}.b,ft_AUC{CVrange_ind}.c),'--','Color',cmp_CVranges(CVrange_ind,:),...
            'HandleVisibility','off')
    end
    xlim([0 max(x_tmp)*1.1])
end

fname = [file_prefix_outputs '/Recruitment_ch' num2str(ch_num_neural)];
saveas(gcf,[fname '.jpg']);
saveas(gcf,[fname '.fig']);

%% Save output data
fname = [file_prefix_outputs '/CAPData_ch' num2str(ch_num_neural) '.mat'];
save(fname,'expdate','stim_rec_dist','ch_num_neural',...
           'trial_values','amp_values_nominal','amp_values_actual',...
           't_window','CV_window','ch_neural_avg','subplot_rows_all','subplot_cols_all',...
           'CAPvsTime_axlim','CAPvsCV_axlim',...
           'CV_bounds','CV_latencybounds','y_CV_bounds','t_bounds_slowest',...
           'AUC','ft_AUC','gof_AUC')
       %'f_window','ch_neural_avg_f','CAPvsFreq_axlim',

%%
%{
%% Plot all recruitment curve data + fits on a single figure
figure
set(gcf,'position',[scrnsz(3)*0.25 scrnsz(4)*0.1 scrnsz(3)*0.6 scrnsz(4)*0.8])
for block_ind = 1:Nblocks
    if block_ind == 1
        start_ind = 1;
    else
        start_ind = Ntrials_block(block_ind-1)+1;
    end
    end_ind = Ntrials_block(block_ind);
    
    x_tmp   = amp_values_actual(start_ind:end_ind);
    
    for CVrange_ind = 1:N_CVranges
        y_tmp   = AUC(CVrange_ind,start_ind:end_ind);
        ft_tmp  = ft_AUC{CVrange_ind,block_ind};
        gof_tmp = ft_AUC{CVrange_ind,block_ind};
        
        subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind); hold on
        plot(x_tmp,y_tmp,'*','color',cmap_PWs(block_ind,:),'DisplayName',[num2str(PW_values_blocks(block_ind)) ' us; r^2 = ' ...
             num2str(gof_AUC{CVrange_ind,block_ind}.rsquare,rsquare_precision)])
        plot(x_values_stimamp,fitfunc_recruitmentcurve(x_values_stimamp,ft_tmp.a,ft_tmp.b,ft_tmp.c),'--', ...
            'color',cmap_PWs(block_ind,:),'HandleVisibility','off')
        
        if gof_AUC{CVrange_ind,block_ind}.rsquare < 0.8
            warning(['Poor recruitment curve fit: ' num2str(PW_values_blocks(block_ind)) ' us (block #' num2str(block_ind) '); '...
                     num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s'])
        end
        
        if block_ind == 1
            if CVrange_ind == 1
                title({expdate,...
                     [num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s']})
            else
                title([num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s'])
            end
            xlabel('Stim amp (mA)')
            ylabel('AUC (\muV*s)')
            legend('location','southeast')
            legend('fontsize',4)
            xlim([0 max(x_tmp)*1.1])
            grid on
        end
    end
end

pause(1) % To get correct axis limits
for CVrange_ind = 1:N_CVranges
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)

    xl = xlim;
    yl = ylim;
    text(xl(2)*0.05,(yl(2)-yl(1))*0.98+yl(1),'PW mask for SD fit, thresh @ 20, 50, 80%: ','fontsize',6)
    text(xl(2)*0.75,(yl(2)-yl(1))*0.98+yl(1),num2str(SDfit_PWmask_FractionalThresh(CVrange_ind,:)),'fontsize',6)
end

for CVrange_ind = 1:N_CVranges
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)

    xl = xlim;
    yl = ylim;
    text(xl(2)*0.05,(yl(2)-yl(1))*0.95+yl(1),'PW mask for max AUC: ','fontsize',6)
    text(xl(2)*0.75,(yl(2)-yl(1))*0.95+yl(1),num2str(maxAUC_PWmask(CVrange_ind,:)),'fontsize',6)
end

for CVrange_ind = 1:N_CVranges
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)

    xl = xlim;
    yl = ylim;
    text(xl(2)*0.05,(yl(2)-yl(1))*0.92+yl(1),'PW mask for SD fit, thresh @ const 20%: ','fontsize',6)
    text(xl(2)*0.75,(yl(2)-yl(1))*0.92+yl(1),num2str(SDfit_PWmask_ConstAUCThresh(CVrange_ind,:)),'fontsize',6)
end

%% Extract thresholds using 20, 50, 80% of max of each recruitment curve 
thresh_targets = [0.2 0.5 0.8];
for block_ind = 1:Nblocks
    for thresh_target_ind = 1:numel(thresh_targets)
        for CVrange_ind = 1:N_CVranges
            y_target_tmp = thresh_targets(thresh_target_ind)*ft_AUC{CVrange_ind,block_ind}.a;
            thresh(CVrange_ind,thresh_target_ind,block_ind) = ...
                fitfunc_recruitmentcurve_inverse(y_target_tmp,ft_AUC{CVrange_ind,block_ind}.a,ft_AUC{CVrange_ind,block_ind}.b,ft_AUC{CVrange_ind,block_ind}.c);

            subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)
            plot(thresh(CVrange_ind,thresh_target_ind,block_ind),y_target_tmp,'^','markersize',15,'color',cmap_PWs(block_ind,:),'HandleVisibility','off')
        end
    end
end

%% Extract thresholds using const AUC thresh for each CV range
thresh_targets_ConstAUC = 0.2;

for thresh_target_ind = 1:numel(thresh_targets_ConstAUC)
    for CVrange_ind = 1:N_CVranges
        % Calculate max AUC and target thresh for each CV range
        y_target_tmp = 0;
        for block_ind = 1:Nblocks
            if maxAUC_PWmask(CVrange_ind,block_ind) == 1
                y_target_tmp = y_target_tmp + ft_AUC{CVrange_ind,block_ind}.a;
            end
        end
        y_target_tmp = y_target_tmp / sum(maxAUC_PWmask(CVrange_ind,:));
        y_target_tmp = y_target_tmp * thresh_targets_ConstAUC(thresh_target_ind);
        
        % Find stim amps at which that threshold is crossed for each CV
        for block_ind = 1:Nblocks
            thresh_constAUC(CVrange_ind,thresh_target_ind,block_ind) = ...
                fitfunc_recruitmentcurve_inverse(y_target_tmp,ft_AUC{CVrange_ind,block_ind}.a,ft_AUC{CVrange_ind,block_ind}.b,ft_AUC{CVrange_ind,block_ind}.c);
            
            subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)
            plot(thresh_constAUC(CVrange_ind,thresh_target_ind,block_ind),y_target_tmp,'s','markersize',15,'color',cmap_PWs(block_ind,:),'HandleVisibility','off')
        end
    end
end

%% Plot and fit Ith-PW
clear h_legendon
Tch_precision = 3;

fo_IthPW_std_log = fitoptions('Method','NonlinearLeastSquares',...
   'Lower',[0.01 0.01],...
   'Upper',[10^6 10^6],...
   'StartPoint',[1 1]);
ft_IthPW_std_log = fittype('log10(rheobase*(1+chronaxie/x))','options',fo_IthPW_std_log); % Standard strength-duration (Lapicque-Weiss) equation

% Fit model to data.
x_values_stimamp = 0:0.1:5000;

figure; set(gcf,'position',[scrnsz(3)*0.25 scrnsz(4)*0.1 scrnsz(3)*0.6 scrnsz(4)*0.8])

% Recruitment curves and SD curves (data + fits) using 20, 50, 80% of max of each recruitment curve 
for thresh_target_ind = 1:numel(thresh_targets)
    for CVrange_ind = 1:N_CVranges
        y_data_tmp = squeeze(thresh(CVrange_ind,thresh_target_ind,:));
        subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind); hold on
        
        if sum(SDfit_PWmask_FractionalThresh(CVrange_ind,:)) > 0
            plot(PW_values_blocks,y_data_tmp,'*-','color',cmp_thresh((thresh_target_ind),:));
            
            [ft_SD{CVrange_ind,thresh_target_ind},gof_SD{CVrange_ind,thresh_target_ind}] = ...
                fit(PW_values_blocks(SDfit_PWmask_FractionalThresh(CVrange_ind,:)==1)',log10(y_data_tmp(SDfit_PWmask_FractionalThresh(CVrange_ind,:)==1)),ft_IthPW_std_log);

            fit_result = ft_SD{CVrange_ind,thresh_target_ind}.rheobase.*(1+ft_SD{CVrange_ind,thresh_target_ind}.chronaxie./x_values_stimamp);

            h_legendon(CVrange_ind,thresh_target_ind) = plot(x_values_stimamp,fit_result,'--','color',cmp_thresh((thresh_target_ind),:),...
                'DisplayName',[num2str(100*thresh_targets(thresh_target_ind)) '%; ' 'r^2 = ' num2str(gof_SD{CVrange_ind,thresh_target_ind}.rsquare,rsquare_precision) ...
                '; Tch = ' num2str(ft_SD{CVrange_ind,thresh_target_ind}.chronaxie,Tch_precision) ' \mus' ...
                '; Irh = ' num2str(ft_SD{CVrange_ind,thresh_target_ind}.rheobase*1000,Tch_precision) ' \muA']);
        else
            h_legendon(CVrange_ind,thresh_target_ind) = plot(PW_values_blocks,y_data_tmp,'*-','color',cmp_thresh((thresh_target_ind),:),...
                'DisplayName',[num2str(100*thresh_targets(thresh_target_ind)) '%']);
        end
        
    end
end

% Recruitment curves and SD data using const AUC thresh for each CV range
for thresh_target_ind = 1:numel(thresh_targets_ConstAUC)
    for CVrange_ind = 1:N_CVranges
        y_data_tmp = squeeze(thresh_constAUC(CVrange_ind,thresh_target_ind,:));
        subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind); hold on
        
        if sum(SDfit_PWmask_ConstAUCThresh(CVrange_ind,:)) > 0
            plot(PW_values_blocks,y_data_tmp,'*-','color','k');
            
            [ft_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind},gof_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}] = ...
                fit(PW_values_blocks(SDfit_PWmask_ConstAUCThresh(CVrange_ind,:)==1)',log10(y_data_tmp(SDfit_PWmask_ConstAUCThresh(CVrange_ind,:)==1)),ft_IthPW_std_log);

            fit_result = ft_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}.rheobase.*(1+ft_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}.chronaxie./x_values_stimamp);

            h_legendon(CVrange_ind,numel(thresh_targets)+thresh_target_ind) = plot(x_values_stimamp,fit_result,'--','color',cmp_thresh((thresh_target_ind),:),...
                'DisplayName',[num2str(100*thresh_targets(thresh_target_ind)) '%; ' 'r^2 = ' num2str(gof_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}.rsquare,rsquare_precision) ...
                '; Tch = ' num2str(ft_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}.chronaxie,Tch_precision) ' \mus' ...
                '; Irh = ' num2str(ft_SD_ConstAUCThresh{CVrange_ind,thresh_target_ind}.rheobase*1000,Tch_precision) ' \muA']);
        else
            h_legendon(CVrange_ind,numel(thresh_targets)+thresh_target_ind) = plot(PW_values_blocks,y_data_tmp,'*-','color','k',...
                'DisplayName',[num2str(100*thresh_targets(thresh_target_ind)) '%']);
        end
    end
end

xlim_min = 0;
xlim_max = max(PW_values_blocks);%x_values_stimamp(end);

for CVrange_ind = 1:N_CVranges
    subplot(subplot_rows_CVranges,subplot_cols_CVranges,CVrange_ind)
    grid on
    if CVrange_ind == 1
        title({expdate,...
             [num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s']})
    else
        title([num2str(CV_bounds(CVrange_ind,1)) ' to ' num2str(CV_bounds(CVrange_ind,2)) ' m/s'])
    end
    xlim([xlim_min xlim_max]); ylim([0 ymax_SD_B])
    xlabel('PW (\mus)')
    ylabel('Ith (mA)')
    legend(h_legendon(CVrange_ind,:));
    legend('fontsize',6)
    
    yl = ylim;
    text(100,yl(2)*0.75,['PW mask for SD fit, thresh @ 20, 50, 80%: ' num2str(SDfit_PWmask_FractionalThresh(CVrange_ind,:))],'fontsize',6)
    text(100,yl(2)*0.7,['PW mask for SD fit, thresh @ const 20%: ' num2str(SDfit_PWmask_ConstAUCThresh(CVrange_ind,:))],'fontsize',6)
end

%% Save data for compiling across experiments
fname = [expdate '_IthPW.mat'];
save(fname,'Nblocks','Ntrials_block','CV_bounds','amp_values_nominal','amp_values_actual',...
     'max_CAP','CV_window','AUC','ft_AUC','gof_AUC',...
     'PW_values','PW_values_blocks','thresh_targets','thresh_targets_ConstAUC','thresh','thresh_constAUC',...
     'ft_SD','gof_SD','ft_SD_ConstAUCThresh','gof_SD_ConstAUCThresh',...
     'SDfit_PWmask_FractionalThresh','maxAUC_PWmask','SDfit_PWmask_ConstAUCThresh')
 
%}