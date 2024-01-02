function atlas_tbl = get_atlas_stats(subject, elec)

run('/home/asm82/vol2surf_Yeo/setup/startup');
fshome = getenv('FREESURFER_HOME');

lh = gifti(['subjects/',subject,'/lh.',elec,'_SCA_fs.func.gii']);
rh = gifti(['subjects/',subject,'/rh.',elec,'_SCA_fs.func.gii']);

lh_emp_map = lh.cdata; rh_emp_map = rh.cdata; 

load('fsaverage_spins.mat');

[~,lh_yeo_labels,ctl] = read_annotation([fshome,'subjects/fsaverage/label/lh.Yeo2011_7Networks_N1000.annot']);
[~,rh_yeo_labels,ctr] = read_annotatiion([fshome,'subjects/fsaverage/label/rh.Yeo2011_7Networks_N1000.annot']);
names = {'VN','SMN','DAN','VAN','LIM','FPN','DMN'};


medial_wall_lh = (ctl.table(1,5) == lh_yeo_labels);
medial_wall_rh = (ctr.table(1,5) == rh_yeo_labels);

lh_emp_map(medial_wall_lh) = NaN;
rh_emp_map(medial_wall_rh) = NaN;

lh_atlas_ind = ctl.table(2:end,5);
rh_atlas_ind = ctr.table(2:end,5);

emp_vals = zeros(1,length(atlas_ind));

for ii = 1:length(atlas_ind)
    atlas_vals = [lh_emp_map(lh_yeo_labels==lh_atlas_ind(ii)); rh_emp_map(rh_yeo_labels==rh_atlas_ind(ii))];
    emp_vals(ii) = median(atlas_vals,'omitnan');
end

num_nulls = length(bigIl(1,:));
null_vals = zeros(num_nulls,length(atlas_ind));

for ii = 1:num_nulls
    for jj = 1:length(atlas_ind)
        lh_null_map = lh_emp_map(bigIl(:,ii));
        rh_null_map = rh_emp_map(bigIr(:,ii));
        null_atlas_vals = [lh_null_map(lh_yeo_labels==lh_atlas_ind(ii)); rh_null_map(rh_yeo_labels==rh_atlas_ind(ii))];
        null_vals(ii,jj) = median(null_atlas_vals,'omitnan');
    end
end

pvals = mean(emp_vals<null_vals);
disp(atlas);
disp(mod);
disp(names);
disp(pvals);

atlas_tbl = array2table([emp_vals;null_vals],'VariableNames',names);
%writetable(atlas_tbl,['atlas_stats_',atlas,'_mod',mod,'.txt']);