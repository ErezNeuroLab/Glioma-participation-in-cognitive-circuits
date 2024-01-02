function get_DAN_tumour_connectivity(subject)

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


% Get tumour fMRI mask
system([fsl,'flirt -in tumour_mask -ref mean_func -init struc2func.mat -applyxfm -o tumour_mask_fmri']);
system([fsl,'fslmaths tumour_mask_fmri -thr 0.5 -bin tumour_mask_fmri']);
system([fsl,'fslmaths tumour_mask_fmri -mas mask_func.nii.gz tumour_mask_fmri']);

% Get DAN fMRI mask
system([fsl,'fslmaths yeo_mprage -thr 1003 -uthr 1003 -bin lh_DAN']);
system([fsl,'fslmaths yeo_mprage -thr 2003 -uthr 2003 -bin rh_DAN']);
system([fsl,'fslmaths lh_DAN -add rh_DAN DAN']);
system([fsl,'flirt -in DAN -ref mean_func -init struc2func.mat -applyxfm -o DAN_fmri_mask']);
system([fsl,'fslmaths DAN_fmri_mask -thr 0.5 -bin DAN_fmri_mask']);

system([fsl, 'fslmaths tumour_mask_fmri -mul DAN_fmri_mask tumour_DAN_overlap']);
system([fsl, 'fslmaths DAN_fmri_mask -sub tumour_DAN_overlap DAN_fmri_mask']);

% Seed-based connectivity of the DAN
system([fsl,'fslmeants -i ffd_clean.nii.gz -m DAN_fmri_mask -o ts.txt']);
T = load('ts.txt');
[img,dims] = read_avw('ffd_clean.nii.gz');
img = reshape(img,dims(1)*dims(2)*dims(3),dims(4));
[out, pval] = corr(T,img');
out(isnan(out)==1) = 0;

out = reshape(out',dims(1),dims(2),dims(3),1);
pval = reshape(pval',dims(1),dims(2),dims(3),1);

% Mask these results to focus on significant results on tumour-infiltrated tissue
tumour_mask = niftiread('tumour_mask_fmri.nii.gz');
num_voxels = sum(sum(sum(tumour_mask)));
thres = 0.05/num_voxels;

out_thres = zeros(size(out));
out_thres(tumour_mask == 1 & pval<thres) = out(tumour_mask == 1 & pval<thres);

% R to Z transform 
out_thres = 0.5*log((1+out_thres)./(1-out_thres));

% Save output image:
INFO = niftiinfo('mean_func.nii.gz');
niftiwrite(single(out_thres),'SCA_maps/DAN_tumour',INFO);
system('gzip -f SCA_maps/*.nii');





