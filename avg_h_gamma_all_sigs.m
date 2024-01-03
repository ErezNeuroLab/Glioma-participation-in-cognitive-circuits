% list of tumor electrodes with significant PSC
list_elecs = {'2EL1','2EL2','2EL3','EL1','EL2'};
list_indices = {5,6,7,5,6};
patient_list = {'2017_06','2017_06','2017_06','2018_05','2018_05'};

average_plots = cell(1,length(patient_list));
average_size = double.empty(0, length(patient_list));
gauss_win = 1000;

for i=1:length(list_elecs)
    elec_plot = get_h_gamma_series(list_elecs{i},list_indices{i},patient_list{i},gauss_win,0);
    average_plots{i} = elec_plot;
    average_size(i) = length(elec_plot);
end
min_val = min(average_size);

% calculate the average for each place between all plots
avg = 0;
for jj=1:length(list_elecs)
    avg = avg + average_plots{jj}(1:min_val);
end

avg = avg/length(list_elecs);
figure
max_time = length(avg)/2000;
time_steps = linspace(0, max_time, length(avg));

% manually calculated the standard error values
ste_vals = arrayfun(@(x)std([average_plots{1}(average_plots{1}==x), average_plots{2}(average_plots{1}==x),average_plots{3}(average_plots{1}==x),average_plots{4}(average_plots{1}==x),average_plots{5}(average_plots{1}==x)]),average_plots{1}(1:min_val))/sqrt(5);
curve1 = avg + ste_vals;
curve2 = avg - ste_vals;
x2 = [time_steps, fliplr(time_steps)];
inBetween = [curve1, fliplr(curve2)];
grayColor = [.8 .8 .8];
ste_vals = fill(x2, inBetween, grayColor);
set(ste_vals,'edgecolor','white');
hold on;
ln = plot(time_steps, avg);
ax = gca;
ax.FontSize = 13;
yline(0,'k--')
xlabel('time (s)', 'FontSize',15,'FontName','Helvetica')
ylabel('% Signal Change', 'FontSize', 15, 'FontName','Helvetica')
ln.LineWidth = 3;
ln.MarkerEdgeColor = 'b';