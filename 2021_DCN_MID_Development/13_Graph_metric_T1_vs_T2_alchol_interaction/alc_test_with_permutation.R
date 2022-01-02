setwd('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov//AUDIT analysis new subid//')
library(lmerTest)
library(dplyr)
# test with fu2 alc group
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


# add audit group
bsl_df=df %>% filter(time==1)
fu2_df=df %>% filter(time==2)
fu2_df=fu2_df %>% mutate(fu2_alc=case_when(audit_total>7~1,
                                          audit_total<8~0))
fu2_df$fu2_alc=as.factor(fu2_df$fu2_alc)
bsl_df$fu2_alc=as.factor(fu2_df$fu2_alc)
df=rbind(bsl_df,fu2_df)


# run net analysis
time_p=list()
time_t=list()
int_t=list()
int_p=list()
var_t=list()
var_p=list()
int_null_t=list()
var_null_t=list()
for (col1 in c('wu_net_aLp','wu_net_adeg','wu_net_aCp')){
  for (col2 in c('fu2_alc')){
    col2show=sprintf('%s-%s',col1,col2)
    tmp_data=df
    f2use=sprintf('%s~mod_pds+handedness+time_f*%s+sex+all_site+age_diff+mean_FD+(1|subcode)',col1,col2)
    lme_model=lmer(as.formula(f2use), data=tmp_data)
    time_t[[col2show]]=summary(lme_model)$coefficients['time_f2','t value']
    time_p[[col2show]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
    var_t[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'t value']
    var_p[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'Pr(>|t|)']
    int_t[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'t value']
    int_p[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'Pr(>|t|)']
    perm_results=node_perm_parallel_function(df,col1,14)
    int_null_t[[col2show]]=perm_results$null_int_t
    var_null_t[[col2show]]=perm_results$null_alc_t
  }
}
net_df=data.frame(time_t=unlist(time_t),
                  time_p=unlist(time_p),
                  var_t=unlist(var_t),
                  var_p=unlist(var_p),
                  int_t=unlist(int_t),
                  int_p=unlist(int_p))
net_df=net_df %>% mutate(fdr_int_p=p.adjust(int_p, method="fdr"))
net_df=net_df %>% mutate(fdr_alc_p=p.adjust(var_p, method="fdr"))
net_df=net_df %>% mutate(fdr_time_p=p.adjust(time_p, method="fdr"))


tmp_quantile=data.frame(l025=unlist(sapply(var_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(var_null_t, function(x){quantile(x,0.975)})))
net_df$alc_l025=tmp_quantile$l025
net_df$alc_u975=tmp_quantile$u975

tmp_quantile=data.frame(l025=unlist(sapply(int_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(int_null_t, function(x){quantile(x,0.975)})))
net_df$int_l025=tmp_quantile$l025
net_df$int_u975=tmp_quantile$u975

net_df=net_df %>% mutate(alc_ifsig=(var_t<alc_l025 | var_t>alc_u975))
net_df=net_df %>% mutate(alc_sig_dir=case_when(alc_ifsig==1 & var_t>0 ~1,
                                               alc_ifsig==1 & var_t<0 ~2,
                                               alc_ifsig==0 ~0))
net_df=net_df %>% mutate(int_ifsig=(int_t<int_l025 | int_t>int_u975))
net_df=net_df %>% mutate(int_sig_dir=case_when(int_ifsig==1 & int_t>0 ~1,
                                               int_ifsig==1 & int_t<0 ~2,
                                               int_ifsig==0 ~0))

net_df=cbind(net_df,net_info)
rownames(net_df)=c('wu_net_aLp','wu_net_adeg','wu_net_aCp')
write.csv(net_df,'alc_net_comp_results.csv')

int_null_t_df=data.frame(int_null_t)
var_null_t=data.frame(var_null_t)
write.csv(int_null_t,'alc_net_time_interaction_null_t.csv')
write.csv(var_null_t,'alc_net_null_t.csv')


plot_model(lme_model,type='int')
summary(lme_model)
d=data.frame(x=df$time,
             y=df$wu_net_aLp,
             id=df$subcode,
             sep=df$fu2_alc)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
longitudinal_plot_by_sep(d,'Time','AUC','')+theme_classic()

# alc vs. control node adeg
time_p=list()
time_t=list()
int_t=list()
int_p=list()
var_t=list()
var_p=list()
int_null_t=list()
var_null_t=list()
for (col1 in colnames(df)[grep('bu_node_adeg_',colnames(df))]){
  for (col2 in c('fu2_alc')){
    col2show=sprintf('%s-%s',col1,col2)
    tmp_data=df
    f2use=sprintf('%s~mod_pds+handedness+time_f*%s+sex+all_site+age_diff+mean_FD+(1|subcode)',col1,col2)
    lme_model=lmer(as.formula(f2use), data=tmp_data)
    time_t[[col2show]]=summary(lme_model)$coefficients['time_f2','t value']
    time_p[[col2show]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
    var_t[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'t value']
    var_p[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'Pr(>|t|)']
    int_t[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'t value']
    int_p[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'Pr(>|t|)']
    perm_results=node_perm_parallel_function(df,col1,14)
    int_null_t[[col2show]]=perm_results$null_int_t
    var_null_t[[col2show]]=perm_results$null_alc_t
  }
}
node_df=data.frame(time_t=unlist(time_t),
                   time_p=unlist(time_p),
                   var_t=unlist(var_t),
                   var_p=unlist(var_p),
                   int_t=unlist(int_t),
                   int_p=unlist(int_p))
node_df=node_df %>% mutate(fdr_int_p=p.adjust(int_p, method="fdr"))
node_df=node_df %>% mutate(fdr_alc_p=p.adjust(var_p, method="fdr"))
node_df=node_df %>% mutate(fdr_time_p=p.adjust(time_p, method="fdr"))


tmp_quantile=data.frame(l025=unlist(sapply(var_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(var_null_t, function(x){quantile(x,0.975)})))
node_df$alc_l025=tmp_quantile$l025
node_df$alc_u975=tmp_quantile$u975

tmp_quantile=data.frame(l025=unlist(sapply(int_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(int_null_t, function(x){quantile(x,0.975)})))
node_df$int_l025=tmp_quantile$l025
node_df$int_u975=tmp_quantile$u975

node_df=node_df %>% mutate(alc_ifsig=(var_t<alc_l025 | var_t>alc_u975))
node_df=node_df %>% mutate(alc_sig_dir=case_when(alc_ifsig==1 & var_t>0 ~1,
                                                 alc_ifsig==1 & var_t<0 ~2,
                                                 alc_ifsig==0 ~0))
node_df=node_df %>% mutate(int_ifsig=(int_t<int_l025 | int_t>int_u975))
node_df=node_df %>% mutate(int_sig_dir=case_when(int_ifsig==1 & int_t>0 ~1,
                                                 int_ifsig==1 & int_t<0 ~2,
                                                 int_ifsig==0 ~0))

node_df=cbind(node_df,node_info)
write.csv(node_df,'alc_node_comp_results.csv')

int_null_t_df=data.frame(int_null_t)
var_null_t=data.frame(var_null_t)
write.csv(int_null_t,'alc_node_time_interaction_null_t.csv')
write.csv(var_null_t,'alc_node_null_t.csv')

node_df %>% filter(int_ifsig==1) %>% select(int_t,int_p,fdr_int_p,int_l025,int_u975,order_idx,estimated_label_full)

node_df %>% filter(alc_ifsig==1) %>% select(var_t,var_p,fdr_alc_p,alc_l025,alc_u975,order_idx,estimated_label_full)







shuffle_alc <- function(cov){
  shuffled_cov=cov
  bsl_alc=cov%>%filter(time==1)%>%select(fu2_alc)
  shuffled_alc=bsl_alc[sample(1:nrow(bsl_alc)),]
  shuffled_cov$shuffled_alc=unlist(list(shuffled_alc,shuffled_alc)) # way to concat two factor
  return(shuffled_cov)
}


node_perm_single_fucntion <- function(gt_df,col){
  gt_df_shuffle_alc=shuffle_alc(gt_df)
  y=gt_df[,col]
  null_model=lmer(y~time_f*shuffled_alc+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df_shuffle_alc)
  null_t=summary(null_model)$coefficients['shuffled_alc1','t value']
  null_int=summary(null_model)$coefficients['time_f2:shuffled_alc1','t value']
  return(cbind(null_t,null_int))
} 
library(doParallel)
node_perm_parallel_function <- function(gt_df,col,ncore2use){
  registerDoParallel(cores=ncore2use)
  tmp_results=foreach (perm=1:1000) %dopar% node_perm_single_fucntion(gt_df,col)
  null_alc_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_t'])))
  null_int_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_int'])))
  result_list=list('null_alc_t'=null_alc_t,
                   'null_int_t'=null_int_t)
  return(result_list)
}

prt=proc.time()
node_perm_parallel_function (gt_df,col,16)
proc.time()-prt
