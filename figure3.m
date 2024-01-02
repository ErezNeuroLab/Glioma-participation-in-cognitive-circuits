% figure 3

run('/imaging/duncan/Ayan/vol2surf_Yeo/setup/startup');
fshome = getenv('FREESURFER_HOME');

[~,lh_yeo_labels,ctl] = read_annotation([fshome,'subjects/fsaverage/label/lh.Yeo2011_7Networks_N1000.annot']);
[~,rh_yeo_labels,ctr] = read_annotation([fshome,'subjects/fsaverage/label/rh.Yeo2011_7Networks_N1000.annot']);

lh_yeo_idx = zeros(size(lh_yeo_labels));rh_yeo_idx = zeros(size(rh_yeo_labels));
lh_labels = ctl.table(:,5);rh_labels = ctr.table(:,5);

count = 0;
for ii = 1:length(lh_labels)
    lh_yeo_idx(lh_yeo_labels==lh_labels(ii)) = count;
    rh_yeo_idx(rh_yeo_labels==rh_labels(ii)) = count;
    count = count+1;
end

chosen_electrodes = {'2017_06', '2EL3', 3, 1; ...
    '2018_04', '2EL3', 3, [5 7]; ...
    '2018_05', 'EL2', [], 4; ...
    '2019_01', '2EL3', 2, [5 7]};

colors_default = 255*ones(7,3); % by default all yeo networks are white

for ii = 1:4
    f1 = vis_fsaverage(chosen_electrodes{ii,1},chosen_electrodes{ii,2});
    
    colors = colors_default;
    red = chosen_electrodes{ii,3};
    blue = chosen_electrodes{ii,4};
    if ~isempty(red)  
        for jj = 1:length(red)
            colors(red,:) = [255, 0, 0]; 
        end
    end
    if ~isempty(blue) 
        for jj = 1:legnth(blue)
            colors(blue,:) = [0, 0, 255];
        end
    end
    
    f2 = CBIG_DrawSurfaceMapsWithBoundary(lh_yeo_idx,rh_yeo_idx,'fsaverage','inflated',min(lh_yeo_idx),max(lh_yeo_idx),colors);
    
    saveas(f1,['plots/',chosen_electrodes{ii,1},'_',chosen_electrodes{ii,2},'_SCA_map,jpg']);
    saveas(f2,['plots/',chosen_electrodes{ii,1},'_',chosen_electrodes{ii,2},'_yeo,jpg']);

end
