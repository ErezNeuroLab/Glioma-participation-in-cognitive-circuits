function stats_all = power_ratio_rotate_stats_tumour(subj_names)
% The function computes stats per electrode based on the rotation approach
% and using instantaneous power using all the data (and not just
% segments/epochs)

% 12/12/2019: I changed the coordinates that are saved in stats_data to the
% shifted (realigned) coordinates

% clear;
% dbstop if error;

%% define paths and load subject data

addpath('/imaging/local/software/fieldtrip/fieldtrip-20160629')
ft_defaults
addpath(genpath('/imaging/duncan/users/ye02/Software/eeglab13_6_5b'))
addpath(genpath('/imaging/duncan/users/ye02/Utils'));

data_folder = '/imaging/duncan/users/ye02/intraop/elecphys_data/';
prm_fn = 'prm.mat';

prm_fn_fp = fullfile(data_folder,prm_fn);
load(prm_fn_fp,'prm'); % load 'prm'

% Get the subject names and electrodes to be included in the analysis.
% [subj_names,elec_names,elec_count] = alt_countF_contrast_inc_electrodes();

% subj_names = {...
%     '2017_02',...
%     '2017_04',...
%     '2017_06',...
%     '2017_07',...
%     '2017_08',...
%     '2017_10',...
%     '2017_11',...
%     '2018_01',...
%     '2018_03',...
%     '2018_04',...
%     '2018_05',...
%     '2018_06',...
%     '2018_07',...
%     '2018_08',...
%     '2019_01',...
%     '2019_02',...
%     };

% specify which rereferenced dataset to analyse
reref_scheme = 'reref_bipolar';
prm.dataset_reref_name = [ '_' reref_scheme]; 

%% Define parameters
prm.spectral_folder = 'spec_analysis'; % Folder for spectral analysis (already defined in prm)
prm.results_folder = '/imaging/duncan/users/ye02/intraop/Results/';
prm.results_power_ratio_stats_folder = ['power_ratio_stats' prm.dataset_reref_name]; %<---- change this line

%sv_dir = fullfile(prm.results_folder,prm.results_power_ratio_stats_folder);
%if ~exist(sv_dir,'dir'); mkdir(sv_dir); end

prm.matlab_all_data_notch50_79_only_fn = ['all_data_notch50_79_only' prm.dataset_reref_name];

prm.power_fn_suffix_notch50_79_only = '_notch50_79_only';
suffix = [prm.power_fn_suffix_notch50_79_only prm.dataset_reref_name];
prm.power_ratio_stats_fn_all_data = ['power_ratio_stats_inst_freq_all_data' suffix];
prm.power_ratio_stats_BrainNet_fn_all_data = ['power_ratio_stats_BrainNet_inst_freq_all_data' suffix];
prm.power_ratio_stats_RenderPlot_fn_all_data = ['power_ratio_stats_render_plot_inst_freq_all_data' suffix];
prm.power_ratio_stats_min_max_vals_fn_all_data = ['power_ratio_stats_min_max_vals_inst_freq_all_data' suffix];
prm.power_ratio_rotate_stats_fn = ['power_ratio_rotate_stats' suffix]; % File name to save for each participant the stats data using the rotation analysis


n_bands = size(prm.bands,1);
bands_names = cell(1,n_bands);
for band_ind = 1:n_bands
    bands_names{band_ind} = sprintf('%d-%d Hz',prm.bands(band_ind,1),prm.bands(band_ind,2));
end


% Contrasts
cont_conds{1} = {'rest','countF'};
cont_conds{2} = {'rest','alt'};
cont_conds{3} = {'countF','alt'};
prm.cont_conds = cont_conds;
n_conts = length(cont_conds);
prm.rotate_n_iterations = 100000;

is_save_perc_sig_change = 1; % 1 - save as % sig change, 0 - save as power ratio

% Contrast data struct
cont_data_struct = struct('cont_conds',[],'data',[],'data_cond_labels',[]);

% is_add_filt_bands_all = [0 1]; % Out of use from 30/8/2019. 1 - add subject-specific bands for filtering from the .m file. 0 - don't
is_add_filt_bands_all = 0; % No additional frequency bands from now on, a notch filter at 79 Hz and harmonics was added for all subjects. 30/8/2019

% Initialize electrode arrays
elec_array_stats_cont_vals = cell(1,n_conts);
elec_array_stats_p_vals = cell(1,n_conts);
for cont_ind = 1 : n_conts
    elec_array_stats_cont_vals{cont_ind} = [];
    elec_array_stats_p_vals{cont_ind} = [];
end % for cond_ind
elec_array_coord = [];
elec_data = [];
tumour_locs = cell(0);

%% Compute power - loop on all subjects
is_use_stats_file = 0; %<--- ?????

for subj_ind = 1 : length(subj_names)
    tic
    fprintf('stats_calc - subject %i\n',subj_ind);
    
    %% Get data
    subj_name = subj_names{subj_ind};
    %this_inc_elec_names = elec_names{subj_ind}; % included electrode names for this subject
    
    % Load matlab file with all data
    matlab_data_folder = fullfile(data_folder,subj_name,prm.matlab_data_folder);
    spectral_folder = fullfile(data_folder,subj_name,prm.spectral_folder);
    data_fn_fp = fullfile(matlab_data_folder,[prm.matlab_all_data_notch50_79_only_fn '_' subj_name]);
    
    % Find tumour electrodes
    tumour = strcmpi({data_fn_fp.data_all.subj_data.elec_data.tumour}, 'Yes');
    ground = (mean(data_fn_fp.data_all.all_rest,2) ==0); 
    this_inc_elec_names = data_fn_fp.data_all.subj_data.elec_names;
    
    % Prepare data for the rotation stats for each of the three
    % contrasts
    cont_data = getContrastDataForStats(data_fn_fp,prm,cont_data_struct);
    
    clear data_all;
    clear stats_data;
    load(data_fn_fp,'data_all');
    
    % stats file name for this subject
    fn_text = [prm.power_ratio_rotate_stats_fn '_' subj_name];
    power_ratio_rotate_stats_fn_fp = fullfile(spectral_folder,fn_text);
    
    if is_use_stats_file == 0
        
        stats_data.subj_data = data_all.subj_data;
        stats_data.cont_conds = prm.cont_conds;
        
        for cont_ind = 1 : n_conts
            
            this_cont = cont_data(cont_ind);
            
            [power_ratio,p_vals] = getRotatePowerStats(...
                this_cont.data,...
                this_cont.data_cond_labels,...
                prm.rotate_n_iterations);
            
            stats_data.power_ratio{cont_ind} = power_ratio;
            stats_data.p_vals{cont_ind} = p_vals;
            
        end % for cont_ind
        
        stats_data.bands = prm.bands;
        stats_data.bands_names = bands_names;
        stats_data.inc_elec_names = this_inc_elec_names;
        % Get electrode coordinates
        nchs = data_all.subj_data.nchs;
        coords = nan(nchs,3);
        for ch_ind = 1 : nchs
            this_elec_data = data_all.subj_data.elec_data(ch_ind);
            coords(ch_ind,:) = [this_elec_data.x_MNI_coord_shift ...
                this_elec_data.y_MNI_coord_shift this_elec_data.z_MNI_coord_shift];
        end
        
        stats_data.elec_coords = coords;
        
        % Save data file for this subject
        %save(power_ratio_rotate_stats_fn_fp,'stats_data','prm');
        
    else
        % Load existing file with stats
        load(power_ratio_rotate_stats_fn_fp,'stats_data','prm');
    end % if is_use_stats_file
    
    % Get relevant data to add to the array of all electrodes across
    % all subjects
    this_subj_tumour_loc = data_all.subj_data.tumour_loc;
    this_elec_tumour_loc = repmat({this_subj_tumour_loc},1,size(stats_data.elec_coords,1));
    
    % Get the indices of the channels that should be included in this
    % analysis for this subject
    inc_chn_inds = getIncChnInds(data_all.subj_data.elec_names,this_inc_elec_names);
    
    % Add to the array of electrodes
    for cont_ind = 1 : n_conts
        elec_array_stats_cont_vals{cont_ind} = cat(1,elec_array_stats_cont_vals{cont_ind},stats_data.power_ratio{cont_ind}(inc_chn_inds,:));
        elec_array_stats_p_vals{cont_ind} = cat(1,elec_array_stats_p_vals{cont_ind},stats_data.p_vals{cont_ind}(inc_chn_inds,:));
    end % for cond_ind
    
    elec_array_coord = cat(1,elec_array_coord,stats_data.elec_coords(inc_chn_inds,:));
    elec_data = cat(2,elec_data,data_all.subj_data.elec_data(inc_chn_inds));
    tumour_locs = cat(2,tumour_locs,this_elec_tumour_loc(inc_chn_inds));
    
    toc
end % for subj_ind

stats_all.cont_conds = prm.cont_conds;
stats_all.elec_array_stats_cont_vals = elec_array_stats_cont_vals;
stats_all.elec_array_stats_p_vals = elec_array_stats_p_vals;
stats_all.elec_array_coord = elec_array_coord;
stats_all.elec_data = elec_data;
stats_all.bands = prm.bands;
stats_all.bands_names = bands_names;
stats_all.rotate_n_iterations = prm.rotate_n_iterations;

% Save power ratio stats for all subjects - only channels that are included in
% the analysis

%fn_fp = fullfile(prm.results_folder,...
%    prm.results_power_ratio_stats_folder,prm.power_ratio_stats_fn_all_data);
%save(fn_fp,'stats_all','prm');

%% Save the updated parameters file
%save(fullfile(prm.data_folder,prm.prm_fn),'prm');

% Get inclusion list, e.g. based on tumour locations
%[inc_elecs_list,inc_elecs_groups_names] = GetElecInclusionList(tumour_locs);

% Save text file with coordinate data and power ratios to be used by
% BrainNet
% SaveForBrainNet(stats_all,prm,inc_elecs_list,inc_elecs_groups_names,...
%     is_save_perc_sig_change,prm.power_ratio_stats_BrainNet_fn_all_data,...
%     prm.power_ratio_stats_min_max_vals_fn_all_data);
%SaveForRenderPlot(stats_all,prm,inc_elecs_list,...
%    inc_elecs_groups_names,prm.power_ratio_stats_RenderPlot_fn_all_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newData = getData(data,exc_start,exc_end,sr)

start_ind = exc_start * sr + 1;
end_ind = size(data,2) - exc_end * sr;

newData = data(:,start_ind:end_ind,:);


%%%%%%%%%%%%%%%%%%%%%%%

function inc_chn_inds = getIncChnInds(elec_names_all,elec_names_inc)
% Get the indices for the channels that should be included in the analysis
% for a given subject

inc_chn_inds = [];
for ind = 1 : length(elec_names_all)
    
    for ind2 = 1 : length(elec_names_inc)
        if strcmp(elec_names_all{ind}, elec_names_inc{ind2})
            inc_chn_inds = [inc_chn_inds ind];
            break;
        end % if strcmp
    end % for ind2
    
end % for ind


%%%%%%%%%%%%%%%%%%%%%%%

function cont_data = getContrastDataForStats(data_fn_fp,prm,cont_data_struct)
% The function prepares the data for the rotation stats for all contrasts
% for a given subject

n_conts = length(prm.cont_conds);

load(data_fn_fp,'data_all');

n_tasks = length(data_all.tasks_order_for_stats);

cont_data = repmat(cont_data_struct,n_conts);

for cont_ind = 1 : n_conts
    this_cont_conds = prm.cont_conds{cont_ind};
    
    cd_data = [];
    cd_data_cond_labels = []; % A vector of size cd_data, with condition label (1,2) for each time point
    for task_ind = 1 : n_tasks
        trial_cond = data_all.tasks_order_for_stats(task_ind).event;
        trial_ind = data_all.tasks_order_for_stats(task_ind).trial;
        if ismember(trial_cond,this_cont_conds)
            switch(trial_cond)
                case 'rest'
                    this_data = getData(data_all.all_rest_bands_power,...
                        prm.epoch_exc_start_end_rest,prm.epoch_exc_start_end_rest,...
                        prm.dn_sample_rate);
                case 'countF'
                    this_data = getData(data_all.all_countF_bands_power{trial_ind},...
                        prm.epoch_exc_start_end_task,prm.epoch_exc_start_end_task,...
                        prm.dn_sample_rate);
                case 'alt'
                    this_data = getData(data_all.all_alt_bands_power{trial_ind},...
                        prm.epoch_exc_start_alt,prm.epoch_exc_start_end_task,...
                        prm.dn_sample_rate);
            end % switch
            
            % Append condition index
            cont_cond_ind = find(strcmp(trial_cond,this_cont_conds)); % condition index for this contrast
            this_cond_labels = ones(1,size(this_data,2)) * cont_cond_ind;
            
            % Append data duration for this condition
            if isempty(cd_data)
                % Set data
                cd_data = this_data;
                cd_data_cond_labels = this_cond_labels;
            else
                % Append data
                cd_data = cat(2,cd_data,this_data);
                cd_data_cond_labels = [cd_data_cond_labels this_cond_labels];
            end
        end % if ismember
    end % for task_ind
    
    % Assign the data to the contrast data
    cont_data(cont_ind).cont_conds = this_cont_conds;
    cont_data(cont_ind).data = cd_data;
    cont_data(cont_ind).data_cond_labels = cd_data_cond_labels;
    
end % for cont_ind

%%%%%%%%%%%%%%%%%%%%%%%

function [power_ratio,p_vals] = getRotatePowerStats(...
    data,data_labels,n_iterations)
% The function computes the stats using the rotation approach for a given
% contrast

nchs = size(data,1);
nfreq = size(data,3);
data_len = size(data,2);

p_vals = nan(nchs,nfreq);

% Compute power ratio for condition2./condition1 in this contrast
data_inds_c1 = data_labels == 1;
data_inds_c2 = data_labels == 2;

data_c1 = (data(:,data_inds_c1,:));
data_c2 = (data(:,data_inds_c2,:));
% Compute power ratio for all channels and frequency bands
power_ratio = squeeze(mean(data_c2,2)) ./ squeeze(mean(data_c1,2));

% Get random data shifts
data_shifts = randi([2 data_len],[1 n_iterations]); % if data_shifts equals 1, the surrogate is equal to the original data, so we exclude that possibility

% Get power ratio estimates for each iteration after shifting the data
power_ratio_surr = nan(nchs,nfreq,n_iterations);

    parfor it_ind = 1 : n_iterations

        this_shift = data_shifts(it_ind);

        % Shift all labels
        this_labels = [data_labels(this_shift:end) data_labels(1:this_shift-1)];

        % Get labels for condition 1 and 2
        this_inds_c1 = this_labels == 1;
        this_inds_c2 = this_labels == 2;

        this_data_c1 = (data(:,this_inds_c1,:));
        this_data_c2 = (data(:,this_inds_c2,:));

        % Compute power ratio for all channels and frequency bands from the
        % rotated data
        power_ratio_surr(:,:,it_ind) = squeeze(mean(this_data_c2,2)) ./ squeeze(mean(this_data_c1,2));

    end % for it_ind


% Compute p values
for ch_ind = 1 : nchs
    for freq_ind = 1 : nfreq
        this_surr = squeeze(power_ratio_surr(ch_ind,freq_ind,:));
        num_larger = length(find(this_surr > power_ratio(ch_ind,freq_ind)));
        p_vals(ch_ind,freq_ind) = num_larger ./ n_iterations;
    end % for freq_ind
end % for ch_ind


%%%%%%%%%%%%%%%%%%%%%%%

function SaveForRenderPlot(stats_all,prm,inc_elecs_list,...
    inc_elecs_groups_names,power_ratio_RenderPlot_fn)


% Get number of electrodes
n_elec = size(stats_all.elec_array_coord,1);

% surface plot structure
sp_str = struct('elec_name',[],'voxels_MNI_x',NaN,...
    'voxels_MNI_y',NaN,'voxels_MNI_z',NaN,'voxels_native_x',NaN,...
    'voxels_native_y',NaN,'voxels_native_z',NaN,...
    'coord_MNI_x',NaN,'coord_MNI_y',NaN,'coord_MNI_z',NaN,...
    'coord_MNI_shift_x',NaN,'coord_MNI_shift_y',NaN,'coord_MNI_shift_z',NaN,...
    'contrast_val',[],'p_val',[],'loc_accuracy',[]);

for cont_ind = 1 : length(prm.cont_conds)
    
    this_cont = prm.cont_conds{cont_ind};
    
    for band_ind = 1 : 1:size(prm.bands,1)
        
        % Get the contrast values for this band
        values = stats_all.elec_array_stats_cont_vals{cont_ind}(:,band_ind); % Get the relevant band data
        
        % Get the percentile values for this band
        percentile_vals = stats_all.elec_array_stats_p_vals{cont_ind}(:,band_ind); % Get the relevant band data
        
        % Reverse the percentiles for decreases in power ratio to get
        % p-values
        inds_small = find(values < 1);
        p_vals = percentile_vals;
        p_vals(inds_small) = 1- p_vals(inds_small);
        
        % Band text for the output file name
        band_txt = [num2str(prm.bands(band_ind,1)) '_' num2str(prm.bands(band_ind,2))]; % Band frequencies text
        
        % Go over all group of electrodes in the inclusion list
        for group_ind = 1 : length(inc_elecs_list)
            
            % Output file
            out_data_file_fp = fullfile(prm.results_folder,...
                prm.results_power_ratio_stats_folder,...
                [power_ratio_RenderPlot_fn '_' this_cont{2} '_' this_cont{1} '_' ...
                band_txt '_' inc_elecs_groups_names{group_ind} '.mat']);
            
            % Go over all electrodes that should be displayed
            counter = 1; % Count number of electrodes for display (for control)
            clear('elecs'); % just in case
            elecs = repmat(sp_str,1,n_elec);
            for elec_ind = 1 : n_elec
                
                this_elec_data = stats_all.elec_data(elec_ind);
                
                % Only display electrodes that should not be excluded
                if isempty(strfind(this_elec_data.exclude,'Yes'))
                    
                    % Apply additional selection criteria on the electrodes
                    if inc_elecs_list{group_ind}(elec_ind) == 1
                        
                        elecs(counter).voxels_MNI_x = this_elec_data.x_MNI_vox_num;
                        elecs(counter).voxels_MNI_y = this_elec_data.y_MNI_vox_num;
                        elecs(counter).voxels_MNI_z = this_elec_data.z_MNI_vox_num;
                        elecs(counter).voxels_native_x = this_elec_data.x_subj_vox_num;
                        elecs(counter).voxels_native_y = this_elec_data.y_subj_vox_num;
                        elecs(counter).voxels_native_z = this_elec_data.z_subj_vox_num;
                        elecs(counter).coord_MNI_x = this_elec_data.x_MNI_coord;
                        elecs(counter).coord_MNI_y = this_elec_data.y_MNI_coord;
                        elecs(counter).coord_MNI_z = this_elec_data.z_MNI_coord;
                        elecs(counter).coord_MNI_shift_x = this_elec_data.x_MNI_coord_shift;
                        elecs(counter).coord_MNI_shift_y = this_elec_data.y_MNI_coord_shift;
                        elecs(counter).coord_MNI_shift_z = this_elec_data.z_MNI_coord_shift;
                        elecs(counter).contrast_val = values(elec_ind);
                        elecs(counter).percentile = percentile_vals(elec_ind);
                        elecs(counter).p_val = p_vals(elec_ind);
                        elecs(counter).elec_name = this_elec_data.elec_name;
                        elecs(counter).loc_accuracy = this_elec_data.localisation_accuracy;
                        
                        % Count number of electrodes for display (for control)
                        counter = counter+1;
                        
                    end % if inc_elecs_list{group_ind}
                    
                end % if isempty(strfind(
            end % for elec_ind
            
            elecs(counter:end) = [];
            
            save(out_data_file_fp,'elecs');
            
        end % for group_ind
        
    end % for band_ind
    
end % for cont_ind