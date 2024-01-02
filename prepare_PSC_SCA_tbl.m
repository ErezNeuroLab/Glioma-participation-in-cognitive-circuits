% prepare R table

load('stats_all.mat');
L = load('stats_2019_01.mat');

patient_id = {stats_all.elec_data.patient_id}';
channel = {stats_all.elec_data.channel}';
cont_countF = stats_all.elec_array_stats_cont_vals{1}(:,6);cont_countF(strcmpi(patient_id,'2019_01'),:) = L.stats_all.elec_array_stats_cont_vals{1}(:,6); cont_countF = (cont_countF-1)*100;
cont_alt = stats_all.elec_array_stats_cont_vals{3}(:,6); cont_alt(strcmpi(patient_id,'2019_01'),:) = nan(9,1); cont_alt = (cont_alt-1)*100;

p_countF = stats_all.elec_array_stats_p_vals{1}(:,6);p_countF(strcmpi(patient_id,'2019_01'),:) = L.stats_all.elec_array_stats_p_vals{1}(:,6); p_countF(p_countF>0.5) = 1 -p_countF(p_countF>0.5); p_countF = p_countF*2;
p_alt = stats_all.elec_array_stats_p_vals{3}(:,6); p_alt(p_alt>0.5) = 1 -p_alt(p_alt>0.5); p_alt = p_alt*2; p_alt(strcmpi(patient_id,'2019_01'),:) = nan(9,1); 

tumour = strcmpi({stats_all.elec_data.tumour},'Yes')';

PSC_tbl = table(patient_id,channel,tumour,cont_countF,cont_alt,p_countF,p_alt);


% SCA_tbl

load('SCA_results');

column_names = {'VN_emp','SMN_emp','DAN_emp','VAN_emp','LIM_emp','FPN_emp','DMN_emp', ...
    'VN_p','SMN_p','DAN_p','VAN_p','LIM_p','FPN_p','DMN_p'};

SCA_tbl = array2table([SCA_results.emp_vals,SCA_results.p_vals],'VariableNames',column_names);

PSC_SCA_tbl = [PSC_tbl,SCA_tbl];

exclude = strcmpi(PSC_SCA_tbl.patient_id,'2018_08');
PSC_SCA_tbl(exclude,:) = [];


% create melted table



