# read AUDIT
bsl_audit=read.csv('/Users/zhipeng/Google Drive/post-doc/IMAGEN_phenotypes/AUDIT/BSL_AUDIT_clean.csv')
fu2_audit=read.csv('/Users/zhipeng/Google Drive/post-doc/IMAGEN_phenotypes/AUDIT/FU2_AUDIT_clean.csv')

# show histogram of AUDIT score
p1=ggplot(bsl_audit,aes(x=audit_total))+
  geom_histogram()+
  xlab('AUDIT')+
  ggtitle('T1 AUDIT total scores (n=2,213)')

p2=ggplot(fu2_audit,aes(x=audit_total))+
  geom_histogram()+
  xlab('AUDIT')+
  ggtitle('T2 AUDIT total scores (n=1,514)')

main=grid.arrange(p1, p2, nrow = 1)
ggsave('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis_drop_out_test/audit_dist.tiff',
       main,width = 8,height = 4)



# drop-out participants
drop_out_all=setdiff(bsl_audit$ID,fu2_audit$ID) 
length(drop_out_all)
bsl012_id=bsl_audit[bsl_audit$audit_total<3,'ID']
drop_out012=intersect(drop_out_all,bsl012_id)
length(bsl012_id)
retain_012=intersect(bsl012_id,fu2_audit$ID)
length(retain_012)

# intersect with T1 available data

# read brain data
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

bsl_brain_df=brain_df %>% filter(time==1)

g_retain_df=bsl_brain_df %>% filter(subcode %in% retain_012) %>% mutate(group_r_rd=1)
dim(g_retain_df)
g_retain_drop_df=bsl_brain_df %>% filter(subcode %in% c(retain_012,drop_out012)) %>% mutate(group_r_rd=2)
dim(g_retain_drop_df)

g_retain_id=g_retain_df %>% select(subcode)
write.csv(g_retain_id,'/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis_drop_out_test/g1_subid_n786.csv')

g_retain_drop_id=g_retain_drop_df %>% select(subcode)
write.csv(g_retain_drop_id,'/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis_drop_out_test/g2_subid_n1054.csv')


# test network properties
var2test_df=rbind(g_retain_df,g_retain_drop_df)
title2show=c('Shortest path length','Strength','Clutering coefficient')
measures=c('wu_net_aLp','wu_net_adeg','wu_net_aCp','behave_behav_large_mean','behave_behav_small_mean','behave_behav_no_mean')
f_vals=list()
p_vals=list()
for (col2test in measures){
tmp=var.test(as.formula(sprintf('%s~group_r_rd',col2test)),data=var2test_df)
f_vals[[col2test]]=as.numeric(tmp$statistic)
p_vals[[col2test]]=as.numeric(tmp$p.value)
}
result_df=data.frame(F_vals=unlist(f_vals),
                     p_vals=unlist(p_vals)) %>% mutate(fdr_p=p.adjust(p_vals,method='fdr'))


write.csv(result_df,'/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis_drop_out_test/drop_out_test_results.csv')

# var(g_retain_df[,'behave_behav_large_mean'])/var(g_retain_drop_df[,'behave_behav_large_mean'])
# 2*(1-pf(0.9607587,df2=785, df1=1053,lower.tail = F))

# test brain response
# variance/sample size/pd




