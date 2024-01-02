% realign_electrodes_nn 

% load coordinates

load('stats_all.mat');
xyz = [stats_all.elec_data.x_subj_vox_num; stats_all.elec_data.y_subj_vox_num; stats_all.elec_data.z_subj_vox_num];
xyz(:,1:6) = [42, 34, 28,36, 33, 32;108, 112, 116, 161, 151, 142; 203, 197, 189, 173, 176, 180];
xyz = correct_elec_coor(xyz,stats_all); % for subject 2018_08 who had a seperate surgery scan
xyz_nn = zeros(size(xyz))';
patient_id = {stats_all.elec_data.patient_id}; % one for each electrode contact

% iterate through each contact, remapping xyz to the closest voxel in the
% cortical reconstruction

for ii = 1:length(patient_id)
    subject_dir = ['subjects/',patient_id{ii}];
    cortex = niftiread([subject_dir, '/yeo_mprage.nii.gz']);
    cortical_voxels = find((cortex>1000 & cortex<1008) | (cortex>2000 & cortex<2008));
    [a,b,c] = ind2sub(size(cortex),cortical_voxels); abc = [a,b,c]; abc = abc-1; % matlab coordinates are off by 1 from voxel coordinates because 0 indices don't exist in matlab
    dist2brain = vecnorm(xyz(:,jj)'-abc,2,2);
    [~,i] = min(dist2brain);
    xyz_nn(ii,:) = abc(i,:);         
end

% update stats_all to include nn'ed coordinates

x_subj_nn = num2cell(xyz_nn(:,1)); y_subj_nn = num2cell(xyz_nn(:,2)); z_subj_nn = num2cell(xyz_nn(:,3));  
[stats_all.elec_data.x_subj_nn] = x_subj_nn{:};
[stats_all.elec_data.y_subj_nn] = y_subj_nn{:};
[stats_all.elec_data.z_subj_nn] = z_subj_nn{:};

save('stats_all.mat','stats_all','prm');
