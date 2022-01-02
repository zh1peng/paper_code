uiopen('F:\Google Drive\post-doc\MID_network_FU2\process_code\activation_LME_new_subid_unbalanced\all_cov_1304_1241.csv',1)

bsl_T=T(T.time==1,:);
fu2_T=T(T.time==2,:);

bsl_subid=pad0(bsl_T.subcode);
fu2_subid=pad0(fu2_T.subcode);

bsl_subcov=bsl_T{:,cov_names};
fu2_subcov=fu2_T{:,cov_names};


bsl_data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_1st_level';
 bsl_out_path='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/prob_maps/BSL_results';
 fu2_data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_1st_level';
 fu2_out_path='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/prob_maps/FU2_results';
save('Prob_maps_BSL_FU2_job_info.mat')
 