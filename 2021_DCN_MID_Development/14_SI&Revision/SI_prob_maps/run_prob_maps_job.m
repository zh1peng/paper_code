addpath(genpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12'));
addpath('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/prob_maps')
load('Prob_maps_BSL_FU2_job_info.mat')
% level_2nd(data_path, output_path,subid, sub_cov,cov_names, con_num)
% level_2nd(bsl_data_path, bsl_out_path,bsl_subid, bsl_subcov,cov_names,2)
% level_2nd(fu2_data_path, fu2_out_path,fu2_subid, fu2_subcov,cov_names,2)

% create_prob_map_function(data_path, output_path,subid, sub_cov,cov_names, con_num,sample_size,rep_time,emask)
create_prob_map_function(bsl_data_path, bsl_out_path,bsl_subid, bsl_subcov,cov_names, 2,50,100,0)
create_prob_map_function(fu2_data_path, fu2_out_path,fu2_subid, fu2_subcov,cov_names, 2,50,100,0)
