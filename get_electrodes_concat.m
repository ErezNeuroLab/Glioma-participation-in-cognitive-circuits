% get maps with all electrodes concatenated

% Set up FSL environment
fslpath = '/usr/local/fsl';
fsl = [fslpath,'/bin/'];

setenv( 'FSLDIR', fslpath);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
path(path,pwd)
clear fsldir fsldirmpath;

%subjects = {'2018_04'}; % CHANGE THIS FOR A NORMAL RUN
subjects = {'2017_06', '2018_04', '2018_05', '2018_08', '2019_01'};
vars = {'countF_sig', 'alt_sig'};
%vars = {'post','shift_nn', 'nn', 'manual'};
%vars = {'nn', 'shift_nn'};
%vars = {'manual'};
%subjects = {'2018_08'};
%'pre', 

% for ii = 1:length(subjects)
% 
%     cd(['subjects/',subjects{ii}]);
%     system([fsl,'fslmerge -t countF_PSC *_countF_PSC.nii.gz']);
%     [~,n] = system([fsl,'fslval countF_PSC dim4']);
%     system([fsl,'fslmaths countF_PSC -Tmean -mul ', num2str(str2num(n)), ' countF_PSC']);
%     
%     system([fsl,'fslmerge -t alt_PSC *_alt_PSC.nii.gz']);
%     [~,n] = system([fsl,'fslval alt_PSC dim4']);
%     system([fsl,'fslmaths alt_PSC -Tmean -mul ', num2str(str2num(n)), ' alt_PSC']);
%     cd('../..');
% end

for jj = 1:length(subjects)
    disp(subjects{jj});
    cd(['subjects/',subjects{jj}]);
    for ii = 1:length(vars)
        
    suffix = vars{ii}
    
    
    %system([fsl,'fslmerge -t pre_shift_spheres *_sphere_highres_pre.nii.gz']);
    %[~,n] = system([fsl,'fslval pre_shift_spheres dim4']);
    %system([fsl,'fslmaths pre_shift_spheres -Tmean -mul ', num2str(str2num(n)), ' pre_shift_spheres']);
    
    %system([fsl,'fslmerge -t post_shift_spheres *_sphere_highres_post.nii.gz']);
    %[~,n] = system([fsl,'fslval post_shift_spheres dim4']);
    %system([fsl,'fslmaths post_shift_spheres -Tmean -mul ', num2str(str2num(n)), ' post_shift_spheres']);
    %cd('../..');
    
    %system([fsl,'fslmerge -t ',suffix,'_shift_spheres *_sphere_highres_',suffix,'.nii.gz']);
    %[~,n] = system([fsl,'fslval ',suffix,'_shift_spheres dim4']);
    %system([fsl,'fslmaths ',suffix,'_shift_spheres -Tmean -mul ', num2str(str2num(n)), ' ',suffix,'_shift_spheres']);
    
    system([fsl,'fslmerge -t ',suffix,' *EL*_',suffix,'.nii.gz']);
    [~,n] = system([fsl,'fslval ',suffix,' dim4']);
    system([fsl,'fslmaths ',suffix,' -Tmean -mul ', num2str(str2num(n)), ' ',suffix]);
    
    
    end
cd('../..');    
    
end
    