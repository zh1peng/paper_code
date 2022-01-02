setwd('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid')
library(lmerTest)
library(dplyr)
library(plyr)
# analysis with gt measures

# merge with all other covs

gt_measures=read.csv('WU_all_gt_measures.csv')
bsl_good_behav=read.csv('bsl_RT_good_1491indx.csv')
fu2_good_behav=read.csv('fu2_RT_good_1365indx.csv')
good_idx_df=rbind(bsl_good_behav,fu2_good_behav)
gt_measures=cbind(gt_measures,good_idx_df) %>% filter(good_behav==1)
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_extraction_new_bu//Nodal analysis//nodes166_info.csv')
bsl_cov=subset(read.csv('bsl_cov1304.csv'),select=-c(X,good_behav))
fu2_cov=subset(read.csv('fu2_cov1240.csv'),select=-c(X,good_behav))
cov_df=rbind(bsl_cov,fu2_cov)
# colnames(cov_df)



gt_df=cbind(cov_df,gt_measures)
numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(gt_df)[grep("node_|net_",colnames(gt_df))])
factor_col=c('subcode','handedness','sex','all_site')

for(col in factor_col){
  gt_df[,col]=as.factor(as.character(gt_df[,col]))
}
for (col in numeric_col){
  gt_df[,col]=as.numeric(as.character(gt_df[,col]))
}
gt_df$time_f=as.factor(gt_df$time) 
#gt_df=subset(gt_df,select=-c(net_adeg,net_abw,net_alocE)) # remove net_adeg as it is keep the same with threshold method


net_p=list()
net_t=list()
net_null_t=list()
# net analysis
for(col in c('net_aCp','net_agE','net_adeg','net_aLp')){
  y=gt_df[,col]
  lme_model=lmer(y~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
  net_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
  net_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
  null_t=numeric()
  for (perm_i in c(1:1000)){
    print(perm_i)
    gt_df_shuffle_time=shuffle_time(gt_df)
    null_model=lmer(y~shuffled_time+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df_shuffle_time)
    null_t=c(null_t,summary(null_model)$coefficients['shuffled_time2','t value'])
  }
  net_null_t[[col]]=null_t
}
net_null_t_df=data.frame(net_null_t)
net_df=data.frame(net_t=unlist(net_t),net_p=unlist(net_p))
net_df=net_df %>% mutate(fdr_p=p.adjust(net_p, method="fdr"))
tmp_quantile=data.frame(l025=unlist(sapply(net_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(net_null_t, function(x){quantile(x,0.975)})))
net_df$l025=tmp_quantile$l025
net_df$u975=tmp_quantile$u975

rownames(net_df)=c('net_aCp','net_agE','net_adeg','net_aLp')

write.csv(net_df,'WU_net_LME.csv')
write.csv(net_null_t_df,'WU_net_LME_null_tvals.csv')

# =======================run permutation=======================
shuffle_time <- function(cov){
  bsl_subcode=as.character(cov$subcode[cov$time==1])
  fu2_subcode=as.character(cov$subcode[cov$time==2])
  # for the subjects that are shown in both BSL and FU2
  # shuffle within subjects to prevent same subject in same wave
  common_subcode=intersect(bsl_subcode,fu2_subcode)
  common_df=cov %>% filter(subcode %in% common_subcode)
  common_df$shuffled_time=ave(common_df$time, common_df$subcode,FUN=sample) 
  # for the subjects that are only shown in BSL or FU2
  # randomly shullfe their time variable
  no_common_subcode=c(setdiff(bsl_subcode,fu2_subcode),setdiff(fu2_subcode,bsl_subcode))
  no_common_df=cov %>% filter(subcode %in% no_common_subcode)
  no_common_df$shuffled_time=sample(no_common_df$time)
  shuffled_cov=rbind(common_df,no_common_df)
  shuffled_cov$shuffled_time=as.factor(shuffled_cov$shuffled_time)
  return(shuffled_cov)
}
write.csv(net_df,'result_net_LME.csv')
write.csv(net_null_t_df,'null_t_values_net_LME.csv')


# note! IT seems like nodal analysis should perform on binary FC instead of weighted FC
# nodal analysis will focus on nodal degree, betweeness
# node_adeg_p=list()
# node_adeg_t=list()
# for(col in colnames(gt_df)[grep('node_adeg_',colnames(gt_df))]){
#   y=gt_df[,col]
#   lme_model=lmer(y~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
#   node_adeg_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
#   node_adeg_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
# }
# 
# node_adeg_df=data_frame(node_adeg_t=unlist(node_adeg_t),node_adeg_p=unlist(node_adeg_p))
# node_adeg_df=node_adeg_df %>% mutate(fdr_p=p.adjust(node_adeg_p, method="fdr"))
# node_adeg_df=cbind(node_adeg_df,node_info)
# write.csv(node_adeg_df,'/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_extraction_new_bu//LME model and plot//node_adeg_LME1.csv')
# 
# node_adeg_df %>% filter(fdr_p<0.05)
# node2test=names(which(node_adeg_df$fdr_p<0.05))
# 
# node_abw_p=list()
# node_abw_t=list()
# for(col in colnames(gt_df)[grep('node_abw_',colnames(gt_df))]){
#   y=gt_df[,col]
#   lme_model=lmer(y~time_f+sex+all_site+handedness+mean_FD+(1|subcode),data=gt_df)
#   node_abw_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
#   node_abw_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
# }
# 
# node_abw_df=data_frame(node_abw_t=unlist(node_abw_t),node_abw_p=unlist(node_abw_p))
# node_abw_df=node_abw_df %>% mutate(fdr_p=p.adjust(node_abw_p, method="fdr"))
# node_abw_df=cbind(node_abw_df,node_info)
# write.csv(node_abw_df,'/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_extraction_new_bu//LME model and plot//node_abw_LME.csv')
# 
# 
# node_abw_df %>% filter(fdr_p<0.05)
# 
# 
# # nodal clustering coefficiency
# node_aCp_p=list()
# node_aCp_t=list()
# for(col in colnames(gt_df)[grep('node_aCp_',colnames(gt_df))]){
#   y=gt_df[,col]
#   lme_model=lmer(y~time_f+sex+all_site+handedness+mean_FD+(1|subcode),data=gt_df)
#   node_aCp_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
#   node_aCp_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
# }
# 
# node_aCp_df=data_frame(node_aCp_t,node_aCp_p)
# node_aCp_df=node_aCp_df %>% mutate(fdr_p=p.adjust(node_aCp_p, method="fdr"))
# node_aCp_df=cbind(node_aCp_df,node_info)
# node_aCp_df %>% filter(fdr_p<0.05)
# 
# # shortest path length
# node_aLp_p=list()
# node_aLp_t=list()
# for(col in colnames(gt_df)[grep('node_aLp_',colnames(gt_df))]){
#   y=gt_df[,col]
#   lme_model=lmer(y~time_f+sex+all_site+handedness+mean_FD+(1|subcode),data=gt_df)
#   node_aLp_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
#   node_aLp_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
# }
# 
# node_aLp_df=data_frame(node_aLp_t,node_aLp_p)
# node_aLp_df=node_aLp_df %>% mutate(fdr_p=p.adjust(node_aLp_p, method="fdr"))
# node_aLp_df=cbind(node_aLp_df,node_info)
# node_aLp_df %>% filter(fdr_p<0.05)







