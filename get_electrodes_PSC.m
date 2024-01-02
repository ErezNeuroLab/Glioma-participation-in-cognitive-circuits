% get electrodes PSC

load('stats_all.mat');
L = load('stats_2019_01.mat');

fsldir = '/usr/local/fsl'; % normally the below code would work without bin
fsl = [fsldir,'/bin/'];
%fsldirmpath = sprintf('%s/etc/matlab',fsldir);
setenv('FSLDIR', fsldir);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');

patient_id = {stats_all.elec_data.patient_id};
channel = {stats_all.elec_data.channel};
cont_countF = stats_all.elec_array_stats_cont_vals{1}(:,6);cont_countF(strcmpi(patient_id,'2019_01'),:) = L.stats_all.elec_array_stats_cont_vals{1}(:,6); cont_countF = (cont_countF-1)*100;
cont_alt = stats_all.elec_array_stats_cont_vals{3}(:,6); cont_alt = (cont_alt-1)*100;

for ii = 1:length(patient_id)
    cd(['subjects/',patient_id{ii}]);
    disp(ii)
    system([fsl,'fslmaths ',channel{ii},'_sphere_highres_shift_nn -mul ', num2str(cont_countF(ii)), ' ', channel{ii},'_countF_PSC_corrected']);
    system([fsl,'fslmaths ',channel{ii},'_sphere_highres_shift_nn -mul ', num2str(cont_alt(ii)), ' ', channel{ii},'_alt_PSC_corrected']);

    %system([fsl,'fslmaths shift_nn_shift_spheres -mul ', num2str(cont_countF(ii)), ' countF_PSC_corrected']);
    %system([fsl,'fslmaths shift_nn_shift_spheres -mul ', num2str(cont_alt(ii)), ' alt_PSC_corrected']);
    cd('../..');
end

subjects = {'2017_06','2018_04','2018_05','2018_08','2019_01'};

for ii = 1:length(subjects)
    cd(['subjects/',subjects{ii}]);
    disp(ii);
    system([fsl,'fslmerge -t countF_PSC_corrected *_countF_PSC_corrected.nii.gz']);
    [~,n] = system([fsl,'fslval countF_PSC_corrected dim4']);
    system([fsl,'fslmaths countF_PSC_corrected -Tmean -mul ', num2str(str2num(n)), ' countF_PSC_corrected']);
    
    system([fsl,'fslmerge -t alt_PSC_corrected *_alt_PSC_corrected.nii.gz']);
    [~,n] = system([fsl,'fslval alt_PSC_corrected dim4']);
    system([fsl,'fslmaths alt_PSC_corrected -Tmean -mul ', num2str(str2num(n)), ' alt_PSC_corrected']);
    cd('../..');
end