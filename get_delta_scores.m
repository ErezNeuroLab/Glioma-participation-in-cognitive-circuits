function [follow_up_scores, subjects] = get_delta_scores(data,names_col)


indices2include = [];

store = find(contains(names_col,'Postop'));
indices2include = [indices2include, store];

postop = names_col(contains(names_col,'Postop'));
postop = postop';
temp = cell2mat(postop);
postop_id = cellstr(temp(:,1:5)); % want just the id


store = find(and(~contains(names_col,postop_id),contains(names_col,'Month3')));
indices2include = [indices2include, store];

month3 = names_col(and(~contains(names_col,postop_id),contains(names_col,'Month3')));
month3 = month3';
temp = cell2mat(month3);
month3_id = cellstr(temp(:,1:5));

log_state = and(and(~contains(names_col,postop_id),~contains(names_col,month3_id)),contains(names_col,'Month12'));


store = find(log_state);
indices2include = [indices2include, store];

month12 = names_col(log_state);
month12 = month12';
temp = cell2mat(month12);
month12_id = cellstr(temp(:,1:5));

follow_up_scores = data(:,indices2include);
subjects = names_col(indices2include);

