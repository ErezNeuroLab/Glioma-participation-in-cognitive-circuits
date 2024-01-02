% get all SCA maps

load('stats_all.mat');

electrodes = {stats_all.elec_data.channel}; %electrodes = electrodes(7:end); % skip 2017_06 for now 
subjects = {stats_all.elec_data.patient_id}; %subjects = subjects(7:end); 

for ii = 1:length(electrodes)
    %if ~exist(['subjects/',subjects{ii},'/SCA_maps/',electrodes{ii},'.nii.gz'], 'file')
        disp(subjects{ii});
        disp(electrodes{ii});
        get_SCA_map_electrode(subjects{ii},electrodes{ii});
    %end
end