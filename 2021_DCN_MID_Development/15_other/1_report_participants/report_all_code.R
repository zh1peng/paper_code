setwd('/Users/zhipeng//Google Drive//post-doc/MID_behavwork_FU2/process_code/activation_LME_new_subid_unbalanced//RT analysis//')

bsl_rt=read.csv('bsl_behav_data1304.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
bsl_cov=read.csv('bsl_cov1304.csv')
fu2_rt=read.csv('fu2_behav_data1241.csv') %>% rename_with( ~ paste("behav", .x, sep = "_"))
fu2_cov=read.csv('fu2_cov1241.csv')
cov=rbind(cbind(bsl_cov,bsl_rt),cbind(fu2_cov,fu2_rt))

# change RT unit to ms
cov$behav_large_mean=cov$behav_large_mean*1000
cov$behav_small_mean=cov$behav_small_mean*1000
cov$behav_no_mean=cov$behav_no_mean*1000
cov$behav_large_minus_no=cov$behav_large_mean-cov$behav_no_mean

cov$behav_left_right_ratio[!is.finite(cov$behav_left_right_ratio)] <- 5.5
library(qwraps2)
library(dplyr)
options(qwraps2_markup = 'markdown')
#options(qwraps2_markup = 'latex')
# cov=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/activation_LME_new_subid_unbalanced//all_cov_1304_1241.csv')
# mean_ci(cov$age)

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
       'BSL PDS'=
         list('Pre-pubertal(1)'  = ~ n_perc0(BL_pds == 1),
              'Beginnig pubertal(2)'  = ~ n_perc0(BL_pds == 2),
              'Mid-pubertal(3)'  = ~ n_perc0(BL_pds == 3),
              'Advanced pubertal(4)'= ~ n_perc0(BL_pds == 4),
              'Post-pubertal(5)'= ~ n_perc0(BL_pds == 5)),
      'Left/Right response'=
        list('mean (sd)'=~mean_sd(behav_left_right_ratio)),
      'RTs'=
        list('Large win'=~mean_sd(behav_large_mean),
             'Small win'=~mean_sd(behav_small_mean),
              'No win'=~mean_sd(behav_no_mean),
             'Large-No win'=~mean_sd(behav_large_minus_no)))
       
#a=summary_table(cov, our_summary1)
sum_T=summary_table(group_by(cov,time),our_summary1)

# add p values

mpvals <-
  sapply(
    list(lm(age ~ time,  data = cov),
         lm(mean_FD ~ time, data = cov)),
    extract_fpvalue)

# Fisher exact test
sex_fpval <- frmtp(fisher.test(table(cov$sex, cov$time))$p.value)
pds_fpval <- frmtp(fisher.test(table(cov$BL_pds, cov$time))$p.value)
handedness_fpval <- frmtp(fisher.test(table(cov$handedness, cov$time))$p.value)

# chi-square
sex_fpval=round(chisq.test(table(cov$sex, cov$time))$p.value,2)
pds_fpval=round(chisq.test(table(cov$BL_pds, cov$time))$p.value,2)
handedness_fpval=round(chisq.test(table(cov$handedness, cov$time))$p.value,2)

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
