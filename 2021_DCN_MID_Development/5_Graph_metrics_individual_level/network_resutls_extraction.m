% extract parameters
function [subs_node, subs_net]=network_resutls_extraction(data_path, subid)
% load('subid_bsl_good_fc_file_1491.mat')
% data_path='/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/fc_analysis/sample_files'
% [bsl_subs_node, bsl_subs_net]=network_resutls_extraction(data_path, subid)
subfiles=fullfile(data_path,strcat('network_sum_bu', subid,'.mat'))
subs_node=node_results_extraction(subfiles);
subs_net=net_results_extraction(subfiles);
end



