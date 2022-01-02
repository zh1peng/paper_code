load('subid_bsl_good_fc_file_1491.mat')
data_path='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_gt_analysis/bsl_fc_code/summary_files'
[bsl_node, bsl_net]=network_resutls_extraction(data_path, subid);

load('subid_fu2_good_fc_file_1365.mat')
data_path='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_gt_analysis/fu2_fc_code/summary_files'
[fu2_node, fu2_net]=network_resutls_extraction(data_path, subid);

clear subid data_path
save('network_results.mat')

