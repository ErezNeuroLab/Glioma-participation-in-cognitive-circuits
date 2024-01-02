function get_SCA_map_electrode(subject, elec)

% Set up FSL environment
fslpath = '/usr/local/fsl';
setenv( 'FSLDIR', fslpath);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
path(path,pwd)
clear fsldir fsldirmpath;

fsl = [fslpath,'/bin/'];

cd(['subjects/',subject]);

if ~exist('yeo_mprage.nii.gz','file')
    system([fsl,'flirt -in Yeo2011_7Networks_N1000 -ref BrainExtractionBrain -applyxfm -usesqform -o yeo_mprage']);
end

% get struc to func transform and mask_func
if ~exist('func2struc.mat','file')
    
system([fsl,'flirt -in BrainExtractionBrain.nii.gz -ref mean_func.nii.gz -omat struc2func.mat -dof 6']);
system([fsl,'fslmaths BrainExtractionBrain.nii.gz -bin brain_mask']);
system([fsl,'flirt -in brain_mask.nii.gz -ref mean_func.nii.gz -init struc2func.mat -applyxfm -o mask_func']);
system([fsl,'fslmaths mask_func.nii.gz -thr 0.5 -bin mask_func.nii.gz']);
system([fsl,'convert_xfm -omat func2struc.mat -inverse struc2func.mat']);

end
% get seed into fMRI space
%system([fsl,'fslmaths ',elec,'_sphere_shift_nn -mas yeo_mprage ',elec,'_sphere_fMRI']); % mask so that ROI is completely within brain
%system([fsl,'fslmaths ',elec,'_sphere_shift_nn -mas brain_mask ',elec,'_sphere_fMRI']);
system([fsl,'fslmaths ',elec,'_sphere_shift_nn_trans -mas brain_mask ',elec,'_sphere_fMRI']); % just for subject 2019_01 electrode 2EL2, which needs to be pushed in more to create an SCA map
system([fsl,'flirt -in ',elec,'_sphere_fMRI -ref mean_func.nii.gz -init struc2func.mat -applyxfm -out ',elec,'_sphere_fMRI']);
system([fsl,'fslmaths ',elec,'_sphere_fMRI -thr 0.5 -bin ',elec,'_sphere_fMRI']);


system([fsl,'fslmeants -i ffd_clean.nii.gz -m ',elec,'_sphere_fMRI.nii.gz -o ts.txt']);
SCA('ts.txt','ffd_clean.nii.gz',elec);

% get it back into structural space
system([fsl,'flirt -in SCA_maps/',elec,' -ref brain_mask -init func2struc.mat -applyxfm -out SCA_maps/',elec]);



cd('../..');