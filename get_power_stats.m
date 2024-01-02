%function get_power_stats(subject)

path2subject = ['/imaging/duncan/users/ye02/intraop/elecphys_data/',subject];
path2data = [path2subject,'/matlab_data/all_data_notch50_79_only_reref_bipolar_',subject,'.mat'];

data = load(path2data);

% concat trials
band = 6;

% countF first
num_trials = length(data.data_all.all_countF_bands_power);
countF = [];

for ii = 1:num_trials
    store = data.data_all.all_countF_bands_power{ii}(:,:,band);
    store = getData(store,data.prm.epoch_exc_start_end_task,data.prm.epoch_exc_start_end_task, data.prm.dn_sample_rate);
    countF = [countF, store];
end
    
% next, alt
num_trials = length(data.data_all.all_alt_bands_power);
alt = [];

for ii = 1:num_trials
    store = data.data_all.all_alt_bands_power{ii}(:,:,band);
    store = getData(store,data.prm.epoch_exc_start_alt,data.prm.epoch_exc_start_end_task, data.prm.dn_sample_rate);
    alt = [alt, store];
end
            
ratio = mean(alt,2) ./ mean(countF,2);


function newData = getData(data,exc_start,exc_end,sr)

start_ind = exc_start * sr + 1;
end_ind = size(data,2) - exc_end * sr;

newData = data(:,start_ind:end_ind,:);
end