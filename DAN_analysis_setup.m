% set up DAN analysis

%list = dir('further_DAN_analysis/20*');

%subjects = {list.name};

patient_data = readtable('patient_data.xlsx');
subjects = patient_data.CBU;

DAN_connectivity = zeros(length(subjects),1);
overlap = zeros(length(subjects),1);

for ii = 1:length(subjects)
    cd(['further_DAN_analysis/',subjects{ii}]);
    map = niftiread('DAN_tumour.nii.gz');
    map = map(map>0);
    
    mask = niftiread('DAN_fmri_mask.nii.gz');
    
    DAN_connectivity(ii) = sum(map) / sum(mask(:));
    
    overlap_map = niftiread('tumour_DAN_overlap.nii.gz');
    overlap(ii) = sum(overlap_map(:));
    
    cd('../..');
end
    
conn_data = table(DAN_connectivity,overlap);
DAN_table = [patient_data, conn_data];
writetable(DAN_table);

    
    