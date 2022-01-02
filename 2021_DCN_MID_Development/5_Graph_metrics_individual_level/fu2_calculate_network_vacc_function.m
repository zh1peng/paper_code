function fu2_calculate_network_vacc_function(subi)
% load subid
% check if the fc file exist
addpath(genpath('/gpfs1/home/z/c/zcao4/matlab_tools/GRETNA-2.0.0_release'))
load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_analysis/fu2_fc_code/subid_fu2_good_fc_file_1365.mat')
data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_1st_level'
fc_file=fullfile(data_path, subid{subi},'fc_map_166.mat');
if exist(fc_file,'file')
tmp=load(fc_file);
M=tmp.con_PPI_final;
M(find(eye(size(M)))) = 0; % chage diognal NaN to 0
M2use={M};
s1=0.05; % search from 5%-40% density
s2=0.4;
deltas=0.02;
n=1000;
% check if the fc file exist
try
[net_sum, node_sum]=gretna_sw_batch_networkanalysis_weight(M2use, s1, s2, deltas, n, 's');
save(['/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_analysis/fu2_fc_code/summary_files/network_sum_',subid{subi},'.mat'],'net_sum','node_sum')
catch
save(['/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_analysis/fu2_fc_code/summary_files/network_sum_',subid{subi},'_error.mat'])
end
end
end

% 