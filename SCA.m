function SCA(ts,D,output)

% Set up FSL environment
setenv( 'FSLDIR', '/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

% Load seed ROI timeseries:
T = load(ts);

% Load BOLD dataset:
[img,dims] = read_avw(D);
img = reshape(img,dims(1)*dims(2)*dims(3),dims(4));

% Perform correlation:
out = corr(T,img');
out = reshape(out',dims(1),dims(2),dims(3),1);
out(isnan(out)==1) = 0;

% Perform r to z transform:
out = 0.5*log((1+out)./(1-out));

% Save output image:
INFO = niftiinfo('mean_func.nii.gz');
niftiwrite(out,['SCA_maps/',output],INFO);

% smooth SCA map
FWHM = 5;
sigma = num2str(FWHM / (2.355));

system('gzip -f SCA_maps/*.nii');
fsl = '/usr/local/fsl/bin/';
%system(['/usr/local/fsl/bin/fslmaths SCA_maps/',output,' -kernel gauss ',sigma,' -fmean SCA_maps/',output]); 
system([fsl,'fslmaths SCA_maps/',output,' -mas mask_func -s ',sigma,' -mas mask_func result1']);
system([fsl,'fslmaths mask_func -s ',sigma,' -mas mask_func result2']);
system([fsl,'fslmaths result1 -div result2 SCA_maps/',output]);
system('gzip -f SCA_maps/*.nii');

%save_avw(out,'SCA_result','f',[2 2 2 1])
%system('/usr/local/fsl/bin/fslcpgeom reg/example_func.nii.gz SCA_result.nii.gz')