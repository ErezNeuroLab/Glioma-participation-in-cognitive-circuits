%function create_elec_spheres(stats_all)

load('stats_all.mat');

fsldir = '/imaging/duncan/users/ye02/fsl'; % normally the below code would work without bin
fsl = [fsldir,'/bin/'];
%fsldirmpath = sprintf('%s/etc/matlab',fsldir);
setenv('FSLDIR', fsldir);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
%path(path, fsldirmpath);
%clear fsldir fsldirmpath;

xyz = [stats_all.elec_data.x_subj_vox_num; stats_all.elec_data.y_subj_vox_num; stats_all.elec_data.z_subj_vox_num];
xyz(:,1:6) = [42, 34, 28,36, 33, 32;108, 112, 116, 161, 151, 142; 203, 197, 189, 173, 176, 180];

RbyC = size(xyz);
patient_id = {stats_all.elec_data.patient_id};
channel = {stats_all.elec_data.channel};
cont_countF = stats_all.elec_array_stats_cont_vals{1}(:,6); cont_countF = (cont_countF-1)*100;
cont_alt = stats_all.elec_array_stats_cont_vals{3}(:,6); cont_alt = (cont_alt-1)*100;
shifts = zeros(RbyC(2),3);
for jj = 1:RbyC(2)
    coord = xyz(:,jj);
    if ~any(isnan(coord))
        cd(['subjects/',patient_id{jj}]);
        
        %x = num2str(xyz(1,jj)); y = num2str(xyz(2,jj)); z =
        %num2str(xyz(3,jj)); sometimes the contact isn't on the brain, so
        %we need to find the voxel on the brain closest to the contacts
        %(which in most cases will be xyz(:,jj)
        
        if strcmpi(patient_id{jj},'2018_08')
            disp(channel{jj});
            x = num2str(xyz(1,jj)); y = num2str(xyz(2,jj)); z = num2str(xyz(3,jj));
            system([fsl,'fslmaths surgery_scan.nii.gz -roi ',x, ' 1 ', y, ' 1 ',z, ' 1 0 1 -bin ',channel{jj},'_surg']);
            system([fsl,'flirt -in ',channel{jj},'_surg -ref BrainExtractionBrain.nii.gz -init surg2preop.mat -applyxfm -o ',channel{jj},' -interp nearestneighbour']);
            elec = niftiread([channel{jj},'.nii.gz']);
            elec_voxel = find(elec);
            [a,b,c] = ind2sub(size(elec),elec_voxel); abc = [a,b,c]; abc = abc-1;
            xyz(:,jj) = abc;
        end
        brain = niftiread('BrainExtractionBrain.nii.gz');
        brain_voxels = find(brain>0);
        %system([fsl,'flirt -in Yeo2011_7Networks_N1000.nii.gz -ref BrainExtractionBrain.nii.gz -applyxfm -usesqform -out yeo_mprage -interp nearestneighbour']);
        %brain = niftiread('Yeo2011_7Networks_N1000.nii.gz');
        
        %brain = niftiread('yeo_mprage.nii.gz');
        %brain = mod(brain,1000);
        %brain_voxels = find((brain>1000 & brain<1008) | (brain>2000 & brain<2008));
        [a,b,c] = ind2sub(size(brain),brain_voxels); abc = [a,b,c]; abc = abc-1; % matlab coordinates are off by 1 from voxel coordinates because 0 indices don't exist in matlab
        dist2brain = vecnorm(xyz(:,jj)'-abc,2,2);
        [~,i] = min(dist2brain);
        shifts(jj,:) = xyz(:,jj)' - abc(i,:);
        
        x = num2str(abc(i,1)); y = num2str(abc(i,2)); z = num2str(abc(i,3));
        
        system([fsl,'fslmaths BrainExtractionBrain.nii.gz -roi ', x, ' 1 ', y, ' 1 ', ...
            z, ' 1 0 1 -bin ', channel{jj}]);
        system([fsl,'fslmaths ',channel{jj},' -kernel sphere 2.5 -fmean -bin ',channel{jj},'_sphere']);
        system([fsl,'flirt -in ',channel{jj},'_sphere -ref ',channel{jj},'_sphere -applyisoxfm 0.5 -out ',channel{jj},'_sphere_highres']);
        system([fsl,'fslmaths ',channel{jj},'_sphere_highres -thr 0.5 -bin ', channel{jj},'_sphere_highres']);
        system([fsl,'fslmaths ',channel{jj},'_sphere_highres -mul ',num2str(cont_countF(jj)),' ',channel{jj},'_countF_PSC']);
        system([fsl,'fslmaths ',channel{jj},'_sphere_highres -mul ',num2str(cont_alt(jj)),' ',channel{jj},'_alt_PSC']);
        
        cd('../..');
    end
end
save('brain_shifts.mat','shifts');
