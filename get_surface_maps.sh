#!/bin/bash

wb=/imaging/local/software/workbench/v1.3.2/bin
FS=/imaging/local/software/freesurfer/6.0.0/x86_64
mesh_path=/imaging/duncan/Ayan/standard_mesh_atlases

export FREESURFER_HOME=$FS
source $FREESURFER_HOME/SetUpFreeSurfer.sh

SUBJECT=$1
ELEC=$2
hemi=$3

cd subjects/${SUBJECT}

$wb/wb_command -surface-apply-affine ${hemi}.pial.surf.gii c_ras.mat ${hemi}.pial_cras.surf.gii

$wb/wb_command -volume-to-surface-mapping SCA_maps/${ELEC}.nii.gz ${hemi}.pial_cras.surf.gii ${hemi}.${ELEC}_SCA.func.gii -trilinear

mris_convert {hemi}.sphere.reg {hemi}.sphere.reg.surf.gii 

$wb/wb_command -metric-resample ${hemi}.${ELEC}_SCA.func.gii ${hemi}.sphere.reg.surf.gii ${mesh_path}/resample_fsaverage/fsaverage_std_sphere.${hemi}.164k_fsavg_${hemi}.surf.gii BARYCENTRIC ${hemi}.${ELEC}_SCA_fs.func.gii

cd ../..

