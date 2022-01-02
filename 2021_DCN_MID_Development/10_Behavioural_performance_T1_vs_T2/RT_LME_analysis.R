library(dplyr)
library(lmerTest)
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

setwd('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/activation_LME_new_subid_unbalanced//RT analysis//')

bsl_rt=read.csv('bsl_behav_data1304.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
bsl_cov=read.csv('bsl_cov1304.csv')
fu2_rt=read.csv('fu2_behav_data1241.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
fu2_cov=read.csv('fu2_cov1241.csv')
df=rbind(cbind(bsl_cov,bsl_rt),cbind(fu2_cov,fu2_rt))
df$behav_large_minus_no=df$behav_large_mean-df$behav_no_mean
df$behav_large_minus_no_in_perc=(df$behav_large_mean-df$behav_no_mean)/df$behav_no_mean*100

numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(df)[grep("behav_",colnames(df))])
factor_col=c('subcode','handedness','sex','all_site')
for(col in factor_col){
  df[,col]=as.factor(as.character(df[,col]))
}
for (col in numeric_col){
  df[,col]=as.numeric(as.character(df[,col]))
}
df$time_f=as.factor(df$time) 

df$behav_left_right_ratio[!is.finite(df$behav_left_right_ratio)] <- 5.5

behav_p=list()
behav_t=list()
behav_null_t=list()
# behav analysis
for(col in c('behav_large_mean','behav_small_mean','behav_no_mean','behav_large_minus_no',  'behav_left_right_ratio')){
#for(col in c('behav_large_minus_no_in_perc')){# add percentage analysis
  y=df[,col]
  lme_model=lmer(y~time_f+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df)
  behav_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
  behav_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
  null_t=numeric()
  for (perm_i in c(1:1000)){
    print(perm_i)
    df_shuffle_time=shuffle_time(df)
    null_model=lmer(y~shuffled_time+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df_shuffle_time)
    null_t=c(null_t,summary(null_model)$coefficients['shuffled_time2','t value'])
  }
  behav_null_t[[col]]=null_t
}
behav_null_t_df=data.frame(behav_null_t)
behav_df=data.frame(behav_t=unlist(behav_t),behav_p=unlist(behav_p))
behav_df=behav_df %>% mutate(fdr_p=p.adjust(behav_p, method="fdr"))
tmp_quantile=data.frame(l025=unlist(sapply(behav_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(behav_null_t, function(x){quantile(x,0.975)})))
behav_df$l025=tmp_quantile$l025
behav_df$u975=tmp_quantile$u975
rownames(behav_df)=c('behav_large_mean','behav_small_mean','behav_no_mean','behav_large_minus_no','behav_left_right_ratio')

write.csv(behav_df,'RT_comparision_results.csv')

d=data.frame(x=df$time,
             y=df$behav_large_mean,
             id=df$subcode,
             sex=df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p1=longitudinal_plot(d,'Time','Large win RT','')+theme_classic()

d=data.frame(x=df$time,
             y=df$behav_small_mean,
             id=df$subcode,
             sex=df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p2=longitudinal_plot(d,'Time','Small win RT','')+theme_classic()

d=data.frame(x=df$time,
             y=df$behav_no_mean,
             id=df$subcode,
             sex=df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p3=longitudinal_plot(d,'Time','No win RT','')+theme_classic()


d=data.frame(x=df$time,
             y=df$behav_large_mean-df$behav_no_mean,
             id=df$subcode,
             sex=df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p4=longitudinal_plot(d,'Time','Large win - No win RT','')+theme_classic()

d=data.frame(x=df$time,
             y=df$behav_left_right_ratio,
             id=df$subcode,
             sex=df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.02)
p5=longitudinal_plot(d,'Time','Left/Right Response Ratio','')+theme_classic()+
  scale_y_continuous(limits=c(0,6)) 


main=grid.arrange(p1,p2,p3,p5,ncol=3)
ggsave('RT_comp_results.tiff',main, width = 8, height = 6)


# add large-no difference analysis
y=df$behav_large_mean-df$behav_no_mean
lme_model=lmer(y~time_f+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df)
behav_null_t=list()
for (perm_i in c(1:1000)){
  print(perm_i)
  df_shuffle_time=shuffle_time(df)
  null_model=lmer(y~shuffled_time+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df_shuffle_time)
  null_t=c(null_t,summary(null_model)$coefficients['shuffled_time2','t value'])
}

# add left/right response ratio analysis 
df_filtered=df %>% filter(behav_left_right_ratio<5)
behav_p=list()
behav_t=list()
behav_null_t=list()
# behav analysis
for(col in c('behav_left_right_ratio')){
  #for(col in c('behav_large_minus_no_in_perc')){# add percentage analysis
  y=df_filtered[,col]
  lme_model=lmer(y~time_f+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df_filtered)
  behav_p[[col]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
  behav_t[[col]]=summary(lme_model)$coefficients['time_f2','t value']
  null_t=numeric()
  for (perm_i in c(1:1000)){
    print(perm_i)
    df_shuffle_time=shuffle_time(df_filtered)
    null_model=lmer(y~shuffled_time+sex+all_site+handedness+mod_pds+age_diff+(1|subcode),data=df_shuffle_time)
    null_t=c(null_t,summary(null_model)$coefficients['shuffled_time2','t value'])
  }
  behav_null_t[[col]]=null_t
}
behav_null_t_df=data.frame(behav_null_t)
behav_df=data.frame(behav_t=unlist(behav_t),behav_p=unlist(behav_p))
behav_df=behav_df %>% mutate(fdr_p=p.adjust(behav_p, method="fdr"))
tmp_quantile=data.frame(l025=unlist(sapply(behav_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(behav_null_t, function(x){quantile(x,0.975)})))
behav_df$l025=tmp_quantile$l025
behav_df$u975=tmp_quantile$u975
rownames(behav_df)=c('behav_left_right_ratio')

write.csv(behav_df,'RT_comparision_results_outlier_removed.csv')


# add large vs. no difference analysis for T1 and T2

df_bsl_cov=df %>% filter(time_f==1) %>% select(sex,all_site,handedness,mod_pds) 
df_bsl_rt1=df %>% filter(time_f==1) %>% select(behav_large_mean) %>% dplyr::rename(rt=behav_large_mean) %>% mutate(rt_g=1)
df_bsl_rt0=df %>% filter(time_f==1) %>% select(behav_no_mean) %>% dplyr::rename(rt=behav_no_mean) %>% mutate(rt_g=0)
df_bsl_analysis=cbind(rbind(df_bsl_cov,df_bsl_cov),rbind(df_bsl_rt0,df_bsl_rt1))
df_bsl_analysis$rt_g=as.factor(df_bsl_analysis$rt_g)
lm_model=lm(rt~rt_g+sex+all_site+handedness+mod_pds,data=df_bsl_analysis)
summary(lm_model)


null_t=numeric()
for (perm_i in c(1:1000)){
  print(perm_i)
  df_bsl_analysis$shuffled_rt_g=sample(df_bsl_analysis$rt_g)
  null_model=lm(rt~shuffled_rt_g+sex+all_site+handedness+mod_pds,data=df_bsl_analysis)
  null_t=c(null_t,summary(null_model)$coefficients['shuffled_rt_g1','t value'])
}
quantile(null_t,0.025)
quantile(null_t,0.975)



df_fu2_cov=df %>% filter(time_f==2) %>% select(sex,all_site,handedness,mod_pds) 
df_fu2_rt1=df %>% filter(time_f==2) %>% select(behav_large_mean) %>% dplyr::rename(rt=behav_large_mean) %>% mutate(rt_g=1)
df_fu2_rt0=df %>% filter(time_f==2) %>% select(behav_no_mean) %>% dplyr::rename(rt=behav_no_mean) %>% mutate(rt_g=0)
df_fu2_analysis=cbind(rbind(df_fu2_cov,df_fu2_cov),rbind(df_fu2_rt0,df_fu2_rt1))
df_fu2_analysis$rt_g=as.factor(df_fu2_analysis$rt_g)
lm_model=lm(rt~rt_g+sex+all_site+handedness+mod_pds,data=df_fu2_analysis)
summary(lm_model)


null_t=numeric()
for (perm_i in c(1:1000)){
  print(perm_i)
  df_fu2_analysis$shuffled_rt_g=sample(df_fu2_analysis$rt_g)
  null_model=lm(rt~shuffled_rt_g+sex+all_site+handedness+mod_pds,data=df_fu2_analysis)
  null_t=c(null_t,summary(null_model)$coefficients['shuffled_rt_g1','t value'])
}
quantile(null_t,0.025)
quantile(null_t,0.975)

