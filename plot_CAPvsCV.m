function plot_CAPvsCV(Ntrials,ch_neural_avg,subplot_rows_all,subplot_cols_all,...
                        CV_window,N_CVranges,CV_bounds,y_CV_bounds,...
                        CV_bounds_slowest,CAPvsCV_axlim,expdate,ch_num_neural,...
                        phase_type,stim_PW,trial_values,amp_values_nominal)

for trial_ind = 1:Ntrials
    % Pull CAP data
    y_tmp = ch_neural_avg(:,trial_ind);

    % Identify the correct subplot for the current trial
    subplot(subplot_rows_all,subplot_cols_all,trial_ind); hold on
    
    % Plot data
    plot(CV_window,y_tmp,'k')
    for CVrange_ind = 1:N_CVranges
        for bound_ind = 1:2 % Lower and upper bounds on each range
            plot(CV_bounds(CVrange_ind,bound_ind).*ones(size(y_CV_bounds)),y_CV_bounds,'k--')
        end
    end
    plot(CV_bounds_slowest.*ones(size(y_CV_bounds)),y_CV_bounds,'r-.')
    axis(CAPvsCV_axlim)
    
    if (trial_ind == 1)
        title({[num2str(expdate) ', ch' num2str(ch_num_neural)],...
               [phase_type ', PW = ' num2str(stim_PW) ' \mus'],...
               ['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA']})
        xlabel('CV [m/s]')
        ylabel('Neural rec [\muV]')
    else
        title(['Trial ' num2str(trial_values(trial_ind)) '; ' ...
               num2str(amp_values_nominal(trial_ind)) ' mA'])
    end
end 