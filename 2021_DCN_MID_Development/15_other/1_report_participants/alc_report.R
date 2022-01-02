df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov//AUDIT analysis new subid//audit_brain_cov.csv')
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/nodes166_info_updated.csv')
numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(df)[grep("bu_node_|net_|wu_node_|behav_|audit_",colnames(df))])
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

bsl_cov=df%>% filter(time==1)
fu2_cov=df%>% filter(time==2)
cov=bsl_cov
cov=fu2_cov


cov$behave_behav_left_right_ratio[!is.finite(cov$behave_behav_left_right_ratio)] <- 5.5
library(qwraps2)
library(dplyr)
options(qwraps2_markup = 'markdown')
our_summary1 <-
  list('Sex' =
         list('Male' = ~n_perc0(sex == 0),
              'Female' =~n_perc0(sex == 1)),
       'Age' =
         list('mean (sd)' = ~ mean_sd(age)),
       'Handedness'=
         list('Left'=~n_perc0(handedness==0),
              'Right'=~n_perc0(handedness==1)),
       'mean FD' =
         list('mean (sd)' = ~ mean_sd(mean_FD)),
       'AUDIT total'=
         list('mean (sd)'=~mean_sd(audit_total)),
       'BSL PDS'=
         list('Pre-pubertal(1)'  = ~ n_perc0(BL_pds == 1),
              'Beginnig pubertal(2)'  = ~ n_perc0(BL_pds == 2),
              'Mid-pubertal(3)'  = ~ n_perc0(BL_pds == 3),
              'Advanced pubertal(4)'= ~ n_perc0(BL_pds == 4),
              'Post-pubertal(5)'= ~ n_perc0(BL_pds == 5)),
       'Left/Right response'=
         list('mean (sd)'=~mean_sd(behave_behav_left_right_ratio)),
       'RTs'=
         list('Large win'=~mean_sd(behave_behav_large_mean),
              'Small win'=~mean_sd(behave_behav_small_mean),
              'No win'=~mean_sd(behave_behav_no_mean)))


sum_T=summary_table(group_by(cov,fu2_alc),our_summary1)


mpvals <-
  sapply(
    list(lm(age ~ fu2_alc,  data = cov),
         lm(mean_FD ~ fu2_alc, data = cov),
         lm(audit_total~fu2_alc, data=cov),
         lm(behave_behav_left_right_ratio~fu2_alc,data=cov)),
    extract_fpvalue)

sex_fpval=chisq.test(table(cov$sex, cov$fu2_alc))$p.value
pds_fpval=chisq.test(table(cov$BL_pds, cov$fu2_alc))$p.value
handedness_fpval=chisq.test(table(cov$handedness, cov$fu2_alc))$p.value


sum_T=cbind(sum_T, "P-value" = "")
sum_T[grepl("mean \\(sd\\)", rownames(sum_T)), "P-value"] <- mpvals
a <- capture.output(print(sum_T))



a[grepl("Sex", a)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(sex_fpval, "|"), a[grepl("Sex", a)])

a[grepl("BSL PDS", a)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(pds_fpval, "|"), a[grepl("BSL PDS", a)])

a[grepl("Handedness", a)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(handedness_fpval, "|"), a[grepl("Handedness", a)])

cat(a, sep = "\n")








# test bsl age
age_pval=wilcox.test(table(cov$age, cov$fu2_alc))
aggregate(age~fu2_alc,data=cov,FUN = 'mean')
ggplot(data=cov,aes(x=fu2_alc,y=age))+
  geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")



wilcox.test(table(cov$age, cov$audit_total))
ggplot(data=cov,aes(x=fu2_alc,y=audit_total))+
  geom_violin()+geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")
