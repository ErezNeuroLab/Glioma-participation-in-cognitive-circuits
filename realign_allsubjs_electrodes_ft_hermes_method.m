% addpath /imaging/ma07/Software/fieldtrip-20190724
% addpath /imaging/local/software/freesurfer/6.0.0/x86_64/matlab
% ft_defaults

%%
clear; clc;
templates_dir =  '/imaging/duncan/users/ye02/intraop/Analysis/Moataz/elec_localization/templates'; %'/imaging/ma07/intraoperative/electrodes_localization/ft_method/templates';
matsv_dir = '/imaging/duncan/users/ye02/intraop/Results/ElecRenderSingleSubj/ft_hermes_method/matlab_data'; %ma07/intraoperative/electrodes_localization/ft_method/saved_data_matlabformat';
figsv_dir = '/imaging/duncan/users/ye02/intraop/Results/ElecRenderSingleSubj/ft_hermes_method/new_realigned_elecs'; %'/imaging/ma07/intraoperative/electrodes_localization/ft_method/figures/gp_avg_elec_render/';
if ~exist(figsv_dir,'dir'); mkdir(figsv_dir); end
if ~exist(matsv_dir,'dir'); mkdir(matsv_dir); end

%% load electrodes coordinates

data_folder = '/imaging/ye02/intraop/Data';
data_fname = 'electrodes_data_MNI';
in_data_file =fullfile(data_folder,data_fname);
load(in_data_file)

%% flip all right hem electrodes to left
all_elecs_MNImm = [ -1*(abs([elec_array.x_MNI_coord]')) [elec_array.y_MNI_coord]' [elec_array.z_MNI_coord]'];

all_elecs_MNImm_noNANs = all_elecs_MNImm;
NAN_rows = isnan(all_elecs_MNImm_noNANs(:,1));
all_elecs_MNImm_noNANs(NAN_rows,:) = [];

%% put electrode coords into fieldtrip structure
load(fullfile(templates_dir, 'SubjectUCI29_elec_acpc_f.mat'));
elecs_ft = elec_acpc_f;
elecs_ft.label =  num2cell(1:size(all_elecs_MNImm_noNANs,1))'; %;{'1','2','3','4'};%{elecs(1:4).elec_name}'; % <--- specify electrode labels here {'1','2','3','4'};
elecs_ft.label = cellfun(@(x) num2str(x), elecs_ft.label,'UniformOutput',false);
elecs_ft.elecpos = all_elecs_MNImm_noNANs; % <--- specify electrode coords here
elecs_ft.cfg.elec.chantype{1} = 'eeg';
elecs_ft.cfg.elec.chanunit{1} = 'V';
elecs_ft.coordsys = 'acpc';

elecs_ft.chanpos = elecs_ft.elecpos;
elecs_ft.tra = eye(length(elecs_ft.elecpos));

elecs_ft.cfg.channel = elecs_ft.label;

elecs_ft.cfg.elec.chanpos = elecs_ft.elecpos;
elecs_ft.cfg.elec.elecpos = elecs_ft.elecpos;
elecs_ft.cfg.elec.chantype = elecs_ft.cfg.elec.chantype(ones(length(elecs_ft.elecpos),1),:);
elecs_ft.cfg.elec.chanunit = elecs_ft.cfg.elec.chanunit(ones(length(elecs_ft.elecpos),1),:);
elecs_ft.cfg.elec.label = elecs_ft.label;

elecs_ft.cfg.elec.cfg.channel = elecs_ft.label;

elecs_ft.cfg.elec.cfg.elec.chanpos = elecs_ft.elecpos;
elecs_ft.cfg.elec.cfg.elec.elecpos = elecs_ft.elecpos;
elecs_ft.cfg.elec.cfg.elec.chantype = elecs_ft.cfg.elec.chantype(ones(length(elecs_ft.elecpos),1),:);
elecs_ft.cfg.elec.cfg.elec.chanunit = elecs_ft.cfg.elec.chanunit(ones(length(elecs_ft.elecpos),1),:);
elecs_ft.cfg.elec.cfg.elec.label = elecs_ft.label;

elecs_ft.cfg.elec.cfg.elec.cfg.channel = elecs_ft.label;
ft_sv_fname = ['elecs_coords_ft_format'];
%     ft_sv_fname = ['elecs_coords_ft_format_' data_fname(25:end)];
%     save(fullfile(matsv_dir, ft_sv_fname),'elecs_ft')

%% get unique number of subjects and contacts for shift correction
elec_names = {elec_array.elec_name};
elec_names_noNANs = elec_names';
elec_names_noNANs(NAN_rows) = [];

for i = 1:size(elec_names_noNANs,1)
    C = strsplit(elec_names_noNANs{i},'_');
    subjs{i,1} = [C{1} '_' C{2}];
    contacts{i,1} = C{3};
end
[unique_subjs,subj_change] = unique(subjs);
[~,~,contact_change] = unique(contacts);
contact_change = find(diff(contact_change) ~=0)+1;

grids = unique([subj_change; contact_change]);

%% correct for brain shifts

elecs_ft_pre_realigned = elecs_ft;
elecs_ft_realigned_hermes = elecs_ft;
elecs_ft_realigned_fsaverage = elecs_ft;

% load(fullfile(templates_dir, 'MNI152_0.8mm_native_hull_lh.mat'))
load(fullfile(templates_dir, 'MNI152_0.8mm_native_hull_cras_lh.mat'))

for g = 1:numel(grids)
    
    elecs_ft_pre_realigned.elecpos = elecs_ft_pre_realigned.elecpos(grids(g):grids(g)+3,:);
    elecs_ft_pre_realigned.chanpos = elecs_ft_pre_realigned.chanpos(grids(g):grids(g)+3,:);
    elecs_ft_pre_realigned.label = elecs_ft_pre_realigned.label(grids(g):grids(g)+3);
    elecs_ft_pre_realigned.tra = eye(length(elecs_ft_pre_realigned.elecpos ));
    
        cfg                         = [];
        cfg.channel            = elecs_ft_pre_realigned.label;
        cfg.keepchannel     = 'yes';
        cfg.elec                  = elecs_ft_pre_realigned;
        cfg.method            = 'headshape';
        cfg.warp                = 'hermes2010';
        cfg.headshape       = hull_lh;
        cfg.feedback          = 'yes';
        
        elecs_ft_realigned_grid_hermes  = ft_electroderealign(cfg);
        elecs_ft_realigned_hermes.elecpos(grids(g):grids(g)+3,:) = elecs_ft_realigned_grid_hermes.elecpos;
        elecs_ft_realigned_hermes.chanpos(grids(g):grids(g)+3,:) = elecs_ft_realigned_grid_hermes.chanpos;
        
        elecs_ft_pre_realigned = elecs_ft;
    
end

%% set up params for fieldtrip surface rendering
face_color = 0.8*[1 1 1];
mesh = [];
[mesh.pos,mesh.tri] = freesurfer_read_surf(fullfile(templates_dir,'MNI_0.8_native_cras_lh.pial'));

f1 = figure(1);clf(1);
subplot(1,2,1)
ft_plot_mesh(mesh,'facecolor',face_color,'facealpha',0.5,'edgealpha',0.5);
ft_plot_sens(elecs_ft_pre_realigned,'facecolor','r','elecshape','sphere');
material dull;
lighting gouraud;
view([-100 20])
camlight;
title('pre shift correc')

subplot(1,2,2)
ft_plot_mesh(mesh,'facecolor',face_color);
ft_plot_sens(elecs_ft_realigned_hermes,'facecolor','r','elecshape','sphere');
material dull;
lighting gouraud;
view([-100 20])
camlight;
title('post shift correc hermes2010')

% print(f1,fullfile(figsv_dir,'gp_avg_plot'),'-r300','-dpng')

%% save realigned electrodes coords in field trip format
ft_sv_realigned_fname = ['ft_format_all_elecs_coords_realigned_hermes'];
elec_names = elec_names_noNANs; %to also save subj names with it in a separate variable
save(fullfile(matsv_dir, ft_sv_realigned_fname),'elecs_ft_realigned_hermes','elec_names','unique_subjs','elecs_ft_pre_realigned')

%% save realigned electrodes coords with subj names
ft_sv_realigned_fname = ['all_elecs_coords_realigned_hermes_w_subjnames'];
k = 0;
for i = 1:length(elec_array)
    if NAN_rows(i) ~= 1
        k = k+1;
        all_elecs_names_shift_correc_coords(i).elecname = elec_array(i).elec_name;
        all_elecs_names_shift_correc_coords(i).x_elecpos = elecs_ft_realigned_hermes.chanpos(k,1);
        all_elecs_names_shift_correc_coords(i).y_elecpos = elecs_ft_realigned_hermes.chanpos(k,2);
        all_elecs_names_shift_correc_coords(i).z_elecpos = elecs_ft_realigned_hermes.chanpos(k,3);
    else
        all_elecs_names_shift_correc_coords(i).elecname = elec_array(i).elec_name;
        all_elecs_names_shift_correc_coords(i).x_elecpos = NaN;
        all_elecs_names_shift_correc_coords(i).y_elecpos = NaN;
        all_elecs_names_shift_correc_coords(i).z_elecpos = NaN;
    end
end

save(fullfile(matsv_dir, ft_sv_realigned_fname),'all_elecs_names_shift_correc_coords')

%% create a surface hull to use for brain shift corrections (already created for MNI template)

% cfg           = [];
% cfg.method    = 'cortexhull';
% cfg.headshape = fullfile(templates_dir,'MNI152_0.8mm.lh.pial');%sv_name;
% cfg.fshome    = '/imaging/ma07/Software/freesurfer_HCPv5.0.6/';
% hull_lh = ft_prepare_mesh(cfg);
% save(fullfile(templates_dir, 'MNI152_0.8mm_native_hull_lh.mat'),'hull_lh')
% 
% % ft_plot_mesh(hull_lh)
% %save hull as .pial
% % load(fullfile(templates_dir,'MNI_0.8_native_hull_lh.mat'))
% % write_surf(fullfile(templates_dir,'MNI_0.8_native_hull_lh.pial'),hull_lh.pos,hull_lh.tri)
% 
% hull_lh = [];
% [hull_lh.pos,hull_lh.tri] = freesurfer_read_surf(fullfile(templates_dir,'MNI_0.8_native_hull_cras_lh.pial'));
% save(fullfile(templates_dir, 'MNI152_0.8mm_native_hull_cras_lh.mat'),'hull_lh')