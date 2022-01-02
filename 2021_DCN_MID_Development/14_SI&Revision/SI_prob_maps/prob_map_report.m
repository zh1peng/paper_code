% Report resutls. Get MNIs for the activiation
% FU -- POS
addpath('F:\Google Drive\zhipeng scripts collection')
mask2use=cellstr(ls('F:\Google Drive\zhipeng scripts collection\fMRI masks\HBM_90+vs+cerebelum+midbrain resliced to spm\*.nii'));
masks2use=fullfile('F:\Google Drive\zhipeng scripts collection\fMRI masks\HBM_90+vs+cerebelum+midbrain resliced to spm',mask2use);
prob_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\FU2_results\all_fwe_pos_binary.nii'
t_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\Normal_2nd_level\FU2_results\spmT_0001.nii';
[FU_pos_results_table,test]=prob_peak_T_in_mask(masks2use,prob_file,t_file);
FU_sig_pos_results=FU_pos_results_table(FU_pos_results_table.above_80>4,:);
% FU -- NEG
prob_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\FU2_results\all_fwe_neg_binary.nii'
t_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\Normal_2nd_level\FU2_results\spmT_0002.nii';
[FU_neg_results_table,test]=prob_peak_T_in_mask(masks2use,prob_file,t_file);
FU_sig_neg_results=FU_neg_results_table(FU_neg_results_table.above_80>4,:);
writetable(FU_sig_pos_results,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\FU2_pos_sig_peaks_5voxel_80.csv','WriteRowNames',true)
writetable(FU_sig_neg_results,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\FU2_neg_sig_peaks_5voxel_80.csv','WriteRowNames',true)



%% BL
prob_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\BSL_results\all_fwe_pos_binary.nii'
t_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\Normal_2nd_level\BSL_results\spmT_0001.nii';
[BL_pos_results_table,test]=prob_peak_T_in_mask(masks2use,prob_file,t_file);
BL_sig_pos_results=BL_pos_results_table(BL_pos_results_table.above_80>4,:);
% FU -- NEG
prob_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\BSL_results\all_fwe_neg_binary.nii'
t_file='F:\Google Drive\post-doc\MID_network_FU2\process_code\Normal_2nd_level\BSL_results\spmT_0002.nii';
[BL_neg_results_table,test]=prob_peak_T_in_mask(masks2use,prob_file,t_file);
BL_sig_neg_results=BL_neg_results_table(BL_neg_results_table.above_80>4,:);
writetable(BL_sig_pos_results,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\BSL_pos_sig_peaks_5voxel_80.csv','WriteRowNames',true)
writetable(BL_sig_neg_results,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\BSL_neg_sig_peaks_5voxel_80.csv','WriteRowNames',true)

%% get all FU data that showed in BL sig
BL_pos_sig_label=BL_sig_pos_results.Properties.RowNames;
FU_pos_all_label=FU_pos_results_table.Properties.RowNames;
[tmp2check, ~,idx]=intersect(BL_pos_sig_label, FU_pos_all_label)
FU_pos2show=FU_pos_results_table(idx,:)
%% check if all FU sig data are in FU_pos2show. In case some regions have more activiation compared with BL
FU_pos_sig_label=FU_sig_pos_results.Properties.RowNames;
FU_pos2show_label=FU_pos_results_table.Properties.RowNames;
num2check=ismember(FU_pos_sig_label,FU_pos2show_label)
sum(num2check)==length(FU_pos_sig_label)
% True!

%% Do the same thing for Negative
BL_neg_sig_label=BL_sig_neg_results.Properties.RowNames;
FU_neg_all_label=FU_neg_results_table.Properties.RowNames;
[tmp2check, ~,idx]=intersect(BL_neg_sig_label, FU_neg_all_label)
FU_neg2show=FU_neg_results_table(idx,:)
FU_neg_sig_label=FU_sig_neg_results.Properties.RowNames;
FU_neg2show_label=FU_neg_results_table.Properties.RowNames;
num2check=ismember(FU_neg_sig_label,FU_neg2show_label)
sum(num2check)==length(FU_neg_sig_label)

writetable(FU_pos2show,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\FU_pos_match2BL_5voxel_80.csv','WriteRowNames',true)
writetable(FU_neg2show,'F:\Google Drive\post-doc\MID_network_FU2\process_code\prob_maps\stats\FU_neg_match2BL_5voxel_80.csv','WriteRowNames',true)
