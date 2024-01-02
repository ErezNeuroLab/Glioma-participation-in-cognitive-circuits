

load('stats_all.mat');

electrodes = {stats_all.elec_data.channel}; %electrodes = electrodes(7:end); % skip 2017_06 for now 
subjects = {stats_all.elec_data.patient_id}; %subjects = subjects(7:end); 

yeo_mean_corr = zeros(length(electrodes), 8);
yeo_network = zeros(length(electrodes),1);

for ii = 1:length(electrodes)
    cd(['subjects/',subjects{ii}]);
    SCA = niftiread(['SCA_maps/',electrodes{ii}, '.nii.gz']);
    %yeo = niftiread('Yeo2011_7Networks_N1000.nii.gz');
    yeo = niftiread('yeo_mprage.nii.gz');
    %yeo = mod(yeo,1000);
    for jj = 1:7
        network_ind = (yeo==(1000+jj)) | (yeo==(2000+jj));
        %network_ind = (yeo==jj);
        yeo_mean_corr(ii,jj) = mean(SCA(network_ind));
    end
    yeo_mean_corr(ii,8) = mean(SCA(yeo>0 & yeo<8));
    
    
    %yeo_mprage = niftiread('yeo_mprage.nii.gz');
    %yeo_mprage = mod(yeo_mprage,1000);
    elec = niftiread([electrodes{ii},'_sphere_shift_nn.nii.gz']);
    %yeo_network(ii) = mode(yeo_mprage(elec>0 & yeo_mprage>0 & yeo_mprage<8));
    yeo_cortex = (yeo>1000 & yeo<1008) | (yeo>2000 & yeo<2008);
    yeo_network(ii) = mode(yeo(elec>0 & yeo_cortex));
    %if yeo_network(ii) == 0
    %    disp('hello');
    %end
    cd('../..');
end

[~,i] = max(yeo_mean_corr');
yeo_connect = i';

