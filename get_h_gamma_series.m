function [PSC_avg_trials] = get_h_gamma_series(elec_num,elec_id, patient_num, gauss_win, istrue)


    load(['/mnt/jane_data/Intraop-Cam/elecphys_data/',patient_num,'/matlab_data/all_data_notch50_79_only_reref_bipolar_',patient_num,'.mat'])
    [hard_series, min_size_hard] = h_gamma_series(data_all.all_alt_bands_power, elec_id);
    [easy_series, min_size_easy] = h_gamma_series(data_all.all_countF_bands_power, elec_id);
    sampling_rate = data_all.sr;

    % smooths the hard trial
    sec_offset = 4;
    start_offset = sampling_rate * sec_offset;
    end_offset = min_size_hard - sampling_rate;
    hard_smooth = smoothdata(hard_series,'gaussian',gauss_win);
    hard_smooth = hard_smooth(start_offset:end_offset);
    
    % smooths the easy trial
    sec_offset = 1;
    start_offset = sampling_rate * sec_offset;
    end_offset = min_size_easy - sampling_rate;
    easy_smooth = smoothdata(easy_series,'gaussian',gauss_win);
    easy_smooth = easy_smooth(start_offset:end_offset);
    
    final_hard = hard_smooth(1:min(length(hard_smooth),length(easy_smooth)));
    final_easy = easy_smooth(1:length(final_hard));
    
    
    % ratio of smoothed hard to easy
    max_time = length(final_easy)/sampling_rate;
    time_steps = linspace(0, max_time, length(final_easy));
    PSC_avg_trials = ((final_hard./final_easy)-1)*100;
    if istrue == 1
        figure
        ln = plot(time_steps, PSC_avg_trials);
        yline(0,'k--')
        xlabel('time (s)')
        ylabel('% Signal Change')
        title(['Ratio of smoothed hard to easy - ', elec_num, ', Patient: ', patient_num])
        ln.LineWidth = 2;
        ln.MarkerEdgeColor = 'b';
    end
end