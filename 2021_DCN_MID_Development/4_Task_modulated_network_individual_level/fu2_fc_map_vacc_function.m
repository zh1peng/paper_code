function []=fu2_fc_map_vacc_function(sub_idx)
tic
addpath('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_map/FU2_code')
addpath(genpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12'))
load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_map/FU2_code/node166_mnis.mat')
load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_map/FU2_code/subid_FU2_1367.mat')
data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_1st_level';
sub_SPM=fullfile(data_path, subid{sub_idx},'SPM.mat');
tmp_SPM=load(sub_SPM);
SPM=tmp_SPM.SPM;
    try
    [~]=mid_fc_map(SPM, mnis,'fc_map_166.mat',sub_idx);
    catch
    error_sub=subid{sub_idx};
    save(['/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_map/BSL_code/error_',num2str(sub_idx),'.mat'], 'error_sub')
    end
    toc
end
   