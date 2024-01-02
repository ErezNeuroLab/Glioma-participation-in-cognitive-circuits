% Set up FSL environment
fslpath = '/usr/local/fsl';
fsl = [fslpath,'/bin/'];

setenv( 'FSLDIR', fslpath);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
path(path,pwd)
clear fsldir fsldirmpath;

% read table 
tbl = readtable('PSC_SCA_tbl.txt');
conts = [tbl.cont_countF, tbl.cont_alt];
p_values = [tbl.p_countF, tbl.p_alt];
conds = {'countF', 'alt'};

for ii = 1:length(tbl.patient_id)
    cd(['subjects/', tbl.patient_id{ii}]);
    for jj = 1:length(conds)
        if (p_values(ii,jj) < 0.05) && (conts(ii,jj) > 0)
            mult = 1;
        elseif (p_values(ii,jj) < 0.05) && (conts(ii,jj) < 0)
            mult = -1;
        else
            mult = 0;
        end
        
        system([fsl,'fslmaths ',tbl.channel{ii},'_sphere_shift_nn -mul ',num2str(mult),' ',tbl.channel{ii},'_', conds{jj}, '_sig']);
    end
    cd ../..
end
        
        