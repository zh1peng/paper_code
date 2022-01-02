library(dplyr)
# merge all the cov

# read BU and WU graph measures
gt_bu_df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/BU_all_gt_measures.csv')
gt_bu_df=gt_bu_df %>% rename_all(function(x) paste0("bu_", x))
gt_wu_df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/WU_all_gt_measures.csv')
gt_wu_df=gt_wu_df %>% rename_all(function(x) paste0("wu_", x))


bsl_good_behav=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/bsl_RT_good_1491indx.csv')
fu2_good_behav=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/fu2_RT_good_1365indx.csv')
good_idx_df=rbind(bsl_good_behav,fu2_good_behav)
gt_measures=cbind(gt_bu_df,gt_wu_df,good_idx_df) %>% filter(good_behav==1) 
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/nodes166_info_updated.csv')


# read ROI activation and fc
bsl_roi=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/bsl_activation_T.csv')
fu2_roi=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/fu2_activation_T.csv')
roi_df=rbind(bsl_roi,fu2_roi) %>% rename_all(function(x) paste0("roi_", x))

bsl_vsl=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/bsl_VS_L_T.csv')
fu2_vsl=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/fu2_VS_L_T.csv')
vsl_df=rbind(bsl_vsl,fu2_vsl) %>% rename_all(function(x) paste0("fc_vsl_", x))

bsl_vsr=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/bsl_VS_R_T.csv')
fu2_vsr=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/activiation_LME//stat and plot//Node166_ROI_extraction/fu2_VS_R_T.csv')
vsr_df=rbind(bsl_vsr,fu2_vsr) %>% rename_all(function(x) paste0("fc_vsr_", x))

# read MGT activation and fc 

bsl_mtg_roi=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/bsl_activation_MTG_T.csv')
fu2_mtg_roi=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/fu2_activation_MTG_T.csv')
mtg_act_df=rbind(bsl_mtg_roi,fu2_mtg_roi) %>% rename_all(function(x) paste0("mgt_roi_", x))

bsl_mtg_vsl=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/bsl_VS_L_MTG_T.csv')
fu2_mtg_vsl=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/fu2_VS_L_MTG_T.csv')
mtg_vsl_df=rbind(bsl_mtg_vsl,fu2_mtg_vsl) %>% rename_all(function(x) paste0("mgt_vsl_", x))


bsl_mtg_vsr=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/bsl_VS_R_MTG_T.csv')
fu2_mtg_vsr=read.csv('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/extract_ROI_with_old_subids/fu2_VS_R_MTG_T.csv')
mtg_vsr_df=rbind(bsl_mtg_vsr,fu2_mtg_vsr) %>% rename_all(function(x) paste0("mgt_vsr_", x))



roi_activiation_fc=cbind(roi_df,vsl_df,vsr_df,mtg_act_df, mtg_vsl_df, mtg_vsr_df, good_idx_df) %>% filter(good_behav==1) %>% select(-c(good_behav))

#read cov
bsl_cov=subset(read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/bsl_cov1304.csv'),select=-c(X,good_behav))
fu2_cov=subset(read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/fu2_cov1240.csv'),select=-c(X,good_behav))
cov_df=rbind(bsl_cov,fu2_cov)

# read rt
bsl_rt=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/activation_LME_new_subid_unbalanced//RT analysis//bsl_behav_data1304.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
fu2_rt=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/activation_LME_new_subid_unbalanced//RT analysis//fu2_behav_data1241.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
rt_df=rbind(bsl_rt,fu2_rt) %>% rename_all(function(x) paste0("behave_", x))







# merge all brain data+ cov
brain_df=cbind(cov_df,roi_activiation_fc,gt_measures,rt_df)
bsl_brain=brain_df%>%filter(time==1)
fu2_brain=brain_df%>%filter(time==2)
common_id=intersect(bsl_brain$subcode,fu2_brain$subcode)



# read audit
audit_cov=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT/final_AUDIT_cov_bsl012.csv')
audit_cov=subset(audit_cov,select = -c(X))


final_common_subcode=intersect(audit_cov$subcode,common_id)

int_audit=audit_cov %>% filter(subcode %in% final_common_subcode) %>% select(starts_with("audit_"))
int_brain=brain_df %>% filter(subcode %in% final_common_subcode)  #%>% rename_all(function(x) paste0("roi_", x))

audit_brain_cov=cbind(int_brain,int_audit)
write.csv(audit_brain_cov,'/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov//AUDIT analysis new subid//audit_brain_cov.csv')




















