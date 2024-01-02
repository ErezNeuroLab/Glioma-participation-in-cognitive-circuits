% pre and post electrode shifts

load('stats_all.mat');

fsldir = '/imaging/duncan/users/ye02/fsl'; % normally the below code would work without bin
fsl = [fsldir,'/bin/'];
%fsldirmpath = sprintf('%s/etc/matlab',fsldir);
setenv('FSLDIR', fsldir);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
%path(path, fsldirmpath);
%clear fsldir fsldirmpath;

xyz_pre_shift = [stats_all.elec_data.x_subj_vox_num; stats_all.elec_data.y_subj_vox_num; stats_all.elec_data.z_subj_vox_num];
xyz_pre_shift(:,1:6) = [42, 34, 28,36, 33, 32;108, 112, 116, 161, 151, 142; 203, 197, 189, 173, 176, 180];

xyz_post_shift = [stats_all.elec_data.x_subj_shift; stats_all.elec_data.y_subj_shift; stats_all.elec_data.z_subj_shift];
xyz_vars = {xyz_pre_shift, '_pre'; xyz_post_shift, '_post'};

patient_id = {stats_all.elec_data.patient_id};
channel = {stats_all.elec_data.channel};

for ii = 1:length(xyz_vars)
    xyz = xyz_vars{ii,1};
    suffix = xyz_vars{ii,2};
for jj = 1:length(channel)
    cd(['subjects/',patient_id{jj}]);
    
    if strcmpi(patient_id{jj},'2018_08')
        disp(channel{jj});
        x = num2str(xyz(1,jj)); y = num2str(xyz(2,jj)); z = num2str(xyz(3,jj));
        system([fsl,'fslmaths surgery_scan.nii.gz -add 1 -bin canvas']);
        system([fsl,'fslmaths canvas -roi ',x, ' 1 ', y, ' 1 ',z, ' 1 0 1 -bin ',channel{jj},'_surg',suffix]);
        system([fsl,'flirt -in ',channel{jj},'_surg -ref BrainExtractionBrain.nii.gz -init surg2preop.mat -applyxfm -o ',channel{jj},suffix,' -interp nearestneighbour']);
        
    end
    
    system([fsl,'fslmaths BrainExtractionBrain.nii.gz -add 1 -bin canvas']);
    system([fsl,'fslmaths canvas -roi ', x, ' 1 ', y, ' 1 ', ...
        z, ' 1 0 1 -bin ', channel{jj},suffix]);
    system([fsl,'fslmaths ',channel{jj},suffix,' -kernel sphere 2.5 -fmean -bin ',channel{jj},'_sphere',suffix]);
    system([fsl,'flirt -in ',channel{jj},'_sphere',suffix,' -ref ',channel{jj},'_sphere',suffix,' -applyisoxfm 0.5 -out ',channel{jj},'_sphere_highres',suffix]);
    system([fsl,'fslmaths ',channel{jj},'_sphere_highres',suffix,' -thr 0.5 -bin ', channel{jj},'_sphere_highres',suffix]);
    %system([fsl,'fslmaths ',channel{jj},'_sphere_highres',suffix,' -mul ',num2str(cont_countF(jj)),' ',channel{jj},'_countF_PSC',suffix]);
    %system([fsl,'fslmaths ',channel{jj},'_sphere_highres',suffix,' -mul ',num2str(cont_alt(jj)),' ',channel{jj},'_alt_PSC',suffix]);
    
    cd('../..');
end
end