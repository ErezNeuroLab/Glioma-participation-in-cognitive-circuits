This folder contains the main code needed to reproduce the findings for
the following manuscript:\
\
Mandal, A.S., Wiener, C., Assem, M., Romero-Garcia, R., Coelho, P.,
McDonald, A., Woodberry, E., Morris, R.C., Price, S.J., Duncan, J.,
Santarius, T., Suckling, J., Hart, M.G., & Erez, Y. (2024).
Tumour-infiltrated cortex participates in large-scale cognitive
circuits. Cortex.\
\
While the raw data cannot be shared due to ethical and privacy
regulations, the code here can be used to conduct similar analyses on
analogous datasets. Any questions can be directed to Ayan Mandal
(Ayan.Mandal@pennmedicine.upenn.edu or email address on personal
website)\
\
[The following scripts were used to determine Percent Signal Change by
contrast:]{.underline}\
\
power_ratio_rotate_stats_tumour.m / get_power stats.m --- uses the
rotation approach to determine stats for each electrode with regards to
their percent signal change for each contrast\
\
[The following scripts were used to visualize the high gamma power
series:]{.underline}\
\
avg_h_gamma_all_sigs.m --- averages the high gamma power series across
subjects\
\
get_h_gama_series.m / h_gamma_series.m --- produces high gamma power
plot for individual subject/electrode\
\
[The following scripts were used to assist with Electrode
Localization:]{.underline}\
\
realign_allsubs_electrodes_ft_hermes_method.m --- projects raw electrode
locations onto the cortical surface\
\
realign_electrodes_nn.m --- maps the electrode locations onto a voxel
coordinate for the functional connectivity analyses\
\
create_elec_spheres.m --- builds spheres around these voxel locations
for both visualization and to be used as seeds in the connectivity
analyses\
\
get_pre_and_post_shifts.m --- another script that does a similar thing
as the one about but more for troubleshooting electrode localization\
\
get_electrodes_PSC.m --- visualizing the electrodes colored by their
Percent Signal change for each contrast\
\
get_sig_elecs.m --- binarizes electrodes as significant vs not
significant PSC for each contrast\
\
get_electrodes_concat.m --- visualizing the electrode locations
themselves, not colored by PSC\
\
[The following scripts were used for the seed-based connectivity
analyses:]{.underline}\
\
SCA.m --- produces seed-based connectivity map give a time-series and
functional connectivity data\
\
get_SCA_map_electrode.m --- calls the above function to produce the
seed-based connectivity map for a given subject/electrode\
\
get_all_SCA_maps.m --- iterates through subjects and electrodes to get
SCA maps using the above function\
\
get_surface_maps.sh --- transforms voxel-based SCA map to surface-based
map\
\
get_atlas_stats.m --- conduct spatial statistical analyses on
surface-based map to determine association with Yeo Networks\
\
get_yeo_stats.m --- not used in final analyses; was used to determine
location of electrode within one of the Yeo Networks\
\
figure3.m --- brings above results together to produce Figure 3 of
manuscript\
\
plot_PSC_SCA.R --- plots relationship between PSC and functional network
connectivity; used for Figure 4\
\
[The following scripts were used for the DAN analyses:]{.underline}\
\
get_DAN_tumour_connectivity,m --- generates DAN connectivity maps of
tumour-infiltrated tissue\
\
DAN_analysis_setup.m --- organizes data for DAN connectivity analyses\
\
get_delta_scores.m --- organizes cognitive data from before and after
surgery\
\
create_DAN_cog_data.m --- organizes DAN and cognitive data together for
further analysis\
\
DAN_location_cog_analyses.R --- conducts statistical analyses on the
relationship between DAN connectivity with tumour location and cognitive
outcomes and also plots these relationships
