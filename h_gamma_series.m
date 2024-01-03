function [full_avg, min_size] = h_gamma_series(location_task, elec_num)

[~,imin]=min(cellfun('length',location_task));
min_size = size(location_task{imin},2);
avg = 0;

for jj = 1:length(location_task)
    avg = avg + (location_task{1,jj}(elec_num,1:min_size,6));
end
full_avg = avg/length(location_task);
end
