library(dplyr)
library(lmerTest)
library(sjPlot)
library(interactions)
library(gridExtra)
library(doParallel)

# ===================================permutation functions=========================================================
shuffle_alc <- function(cov){
  shuffled_cov=cov
  bsl_alc=cov%>%filter(time==1)%>%select(alc_changes)
  shuffled_alc_changes=bsl_alc[sample(1:nrow(bsl_alc)),]
  shuffled_cov$shuffled_alc_changes=unlist(list(shuffled_alc_changes,shuffled_alc_changes)) # way to concat two factor
  return(shuffled_cov)
}

perm_alc_changes_single_fucntion <- function(gt_df,col){
  gt_df_shuffle_alc=shuffle_alc(gt_df)
  y=gt_df[,col]
  null_model=lmer(y~time_f*shuffled_alc_changes+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df_shuffle_alc)
  null_t=summary(null_model)$coefficients['shuffled_alc_changes','t value']
  null_int=summary(null_model)$coefficients['time_f2:shuffled_alc_changes','t value']
  return(cbind(null_t,null_int))
} 

perm_alc_changes_parallel_function <- function(gt_df,col,ncore2use){
  registerDoParallel(cores=ncore2use)
  tmp_results=foreach (perm=1:1000) %dopar% perm_alc_changes_single_fucntion(gt_df,col)
  null_alc_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_t'])))
  null_int_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_int'])))
  result_list=list('null_alc_t'=null_alc_t,
                   'null_int_t'=null_int_t)
  return(result_list)
}
# ==============================================================================================================




# read and do formating
df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov//AUDIT analysis new subid//audit_brain_cov.csv')
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/nodes166_info_updated.csv')
numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(df)[grep("bu_|wu_|behav_|audit_",colnames(df))])
factor_col=c('subcode','handedness','sex','all_site')

for(col in factor_col){
  df[,col]=as.factor(as.character(df[,col]))
}
for (col in numeric_col){
  df[,col]=as.numeric(as.character(df[,col]))
}
df$time_f=as.factor(df$time) 

# add bsl audit
df$audit_total_bsl=c(df[df$time==1,"audit_total"],df[df$time==1,"audit_total"])

#df=df %>% filter(audit_total_bsl<2)

# add audit group
bsl_df=df %>% filter(time==1)
fu2_df=df %>% filter(time==2)

# calculate alc changes
alc_changes=fu2_df$audit_total-bsl_df$audit_total
bsl_df$alc_changes=as.numeric(alc_changes)
fu2_df$alc_changes=as.numeric(alc_changes)


df=rbind(bsl_df,fu2_df)



# run net+mtg_vsr analysis
time_p=list()
time_t=list()
int_t=list()
int_p=list()
var_t=list()
var_p=list()
int_null_t=list()
var_null_t=list()
for (col1 in c('wu_net_aLp','wu_net_adeg','wu_net_aCp','mgt_vsr_node_MTG')){
  for (col2 in c('alc_changes')){
    col2show=sprintf('%s-%s',col1,col2)
    tmp_data=df
    f2use=sprintf('%s~mod_pds+handedness+time_f*%s+sex+all_site+age_diff+mean_FD+(1|subcode)',col1,col2)
    lme_model=lmer(as.formula(f2use), data=tmp_data)
    time_t[[col2show]]=summary(lme_model)$coefficients['time_f2','t value']
    time_p[[col2show]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
    var_t[[col2show]]=summary(lme_model)$coefficients[sprintf('%s',col2),'t value']
    var_p[[col2show]]=summary(lme_model)$coefficients[sprintf('%s',col2),'Pr(>|t|)']
    int_t[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s',col2),'t value']
    int_p[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s',col2),'Pr(>|t|)']
    perm_results=perm_alc_changes_parallel_function(df,col1,14)
    int_null_t[[col2show]]=perm_results$null_int_t
    var_null_t[[col2show]]=perm_results$null_alc_t
  }
}
result_df=data.frame(time_t=unlist(time_t),
                     time_p=unlist(time_p),
                     var_t=unlist(var_t),
                     var_p=unlist(var_p),
                     int_t=unlist(int_t),
                     int_p=unlist(int_p))




tmp_quantile=data.frame(l025=unlist(sapply(var_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(var_null_t, function(x){quantile(x,0.975)})))
result_df$alc_l025=tmp_quantile$l025
result_df$alc_u975=tmp_quantile$u975

tmp_quantile=data.frame(l025=unlist(sapply(int_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(int_null_t, function(x){quantile(x,0.975)})))
result_df$int_l025=tmp_quantile$l025
result_df$int_u975=tmp_quantile$u975

result_df=result_df %>% mutate(alc_ifsig=(var_t<alc_l025 | var_t>alc_u975))
result_df=result_df %>% mutate(alc_sig_dir=case_when(alc_ifsig==1 & var_t>0 ~1,
                                                     alc_ifsig==1 & var_t<0 ~2,
                                                     alc_ifsig==0 ~0))
result_df=result_df %>% mutate(int_ifsig=(int_t<int_l025 | int_t>int_u975))
result_df=result_df %>% mutate(int_sig_dir=case_when(int_ifsig==1 & int_t>0 ~1,
                                                     int_ifsig==1 & int_t<0 ~2,
                                                     int_ifsig==0 ~0))


write.csv(result_df,'/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis/alc_changes_results.csv')



# make plot
# p1 audit plot
# p2 vs ppi interaction
# p3 network interaction

p1=ggplot(bsl_df, aes(x=audit_changes))+geom_histogram(binwidth=1)+
  xlab('AUDIT changes')+
  # geom_segment(aes(x=0.59,xend=0.59,y=0,yend=100),color='coral1')+
  # geom_segment(aes(x=4.83,xend=4.83,y=0,yend=100),color='limegreen')+
  # geom_segment(aes(x=9.08,xend=9.08,y=0,yend=100),color='cornflowerblue')+
  theme_classic()



lme_model=lmer(mgt_vsr_node_MTG~time_f*alc_changes+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
p2=plot_model(lme_model,type='int',mdrt.values = 'meansd')+
  ylab('Marginal Effects')+
  xlab('Time')+
  ggtitle('Left middle temporal gyrus')+
  scale_colour_discrete(name='AUDIT changes')+
  theme_classic()+theme(legend.position = "none")


lme_model=lmer(wu_net_aLp~time_f*alc_changes+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
p3=plot_model(lme_model,type='int',mdrt.values = 'meansd')+
  ylab('Marginal Effects')+
  xlab('Time')+
  ggtitle('Shortest path length')+
  scale_colour_discrete(name='AUDIT changes')+
  theme_classic()+
  theme(legend.position = c(0.8,0.75),
        legend.background = element_rect(color='transparent',fill=NA),
        legend.key=element_rect(color='transparent',fill=NA))

main=grid.arrange(p1, p2,p3, nrow = 1)
ggsave('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis/alc_changes_interaction_plot.tiff',
       main,width = 10,height = 4)