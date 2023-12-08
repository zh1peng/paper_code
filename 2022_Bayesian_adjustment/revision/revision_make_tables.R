
library(lmerTest)
library(purrr)
library(broom)
library(tidyr)
library(dplyr)

cbind.fill <- function(...){ # can use list instead
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

setwd('F:/Google Drive/post-doc/Bayesian_Project/new_model/revision/gamma_narrower')
filenames=c('L_ct', 'R_ct')
all_site_es=read.csv('all_site_es.csv')
all_site_sd=read.csv('all_site_sd.csv')

all_ROImean=data.frame()
all_ROImode=data.frame()
all_ROIrhat=data.frame()
all_ROIess=data.frame()
all_ROIlower=data.frame()
all_ROIupper=data.frame()

for (file_i in filenames){
  load(sprintf('Posterior_extraction_%s.RData',file_i))
  all_ROImean=cbind.fill(all_ROImean,all_mean)
  all_ROImode=cbind.fill(all_ROImode,all_mode)
  all_ROIrhat=cbind.fill(all_ROIrhat,all_rhat)
  all_ROIess=cbind.fill(all_ROIess,all_ess)
  all_ROIlower=cbind.fill(all_ROIlower,all_lower)
  all_ROIupper=cbind.fill(all_ROIupper,all_upper)
}

#############################################################
# scatter plot for mode value vs. raw value (shrinkage plot)
# using mean/median seems better than mode for this plot
#############################################################
site_n=as.numeric(readRDS('F:/Google Drive/post-doc/Bayesian_Project/new_model/site_N.rds')) #
site_n_pad=c(9999, site_n) # add a large site_n for mu
col2test=colnames(all_site_es)[grep('.*\\_thickavg$',colnames(all_site_es))]
flatten_bayes_es_mean=c(all_ROImean[3:23,col2test]) # mean
flatten_bayes_es_mode=c(all_ROImode[3:23,col2test]) # mode
flatten_bayes_overes_upper=c(rep(all_ROIupper[2,col2test],each=21))
flatten_bayes_overes_lower=c(rep(all_ROIlower[2,col2test],each=21))
flatten_es=c(as.matrix(all_site_es[,col2test]))
matrix_n=matrix(rep(site_n,68),ncol=68)
flatten_n=c(matrix_n)
flatten_sd=c(as.matrix(all_site_sd[,col2test]))
sorted_n=sort(site_n,index.return=T) 

study_id=rep(paste0('study_',c(1:21))[order(sorted_n$ix)],68) # make orders match study1 -> n=27
roi=rep(col2test,each=21)
df=data.frame(orig_es=flatten_es,
              bayes_mean_es=flatten_bayes_es_mean,
              bayes_mode_es=flatten_bayes_es_mode,
              bayes_overes_upper=flatten_bayes_overes_upper,
              bayes_overes_lower=flatten_bayes_overes_lower,
              sd=flatten_sd,
              sample_size=flatten_n,
              roi=roi,
              study_id=study_id)%>%
              mutate(mode_adj=abs(orig_es-bayes_mode_es),
                     mean_adj=abs(orig_es-bayes_mean_es),
                     overes_hdi=bayes_overes_upper-bayes_overes_lower) 


# test if the adjustment greater than 0
study_n_df=data.frame(study_id=paste0('study_',c(1:21))[order(sorted_n$ix)],
                      sample_n=site_n)
res_df=df %>% 
  select(orig_es,bayes_mode_es,sample_size,roi,study_id) %>%
  mutate(dist_es=abs(orig_es-bayes_mode_es)) %>% 
  group_by(study_id) %>% 
  nest() %>% 
  mutate(t_model=map(data,~t.test(.x$dist_es,mu=0,alternative='greater'))) %>%  #
  mutate(tidy_db =  map(t_model,tidy)) %>% # 
  unnest(tidy_db) %>% 
  select(study_id,estimate, statistic, p.value, parameter) %>% 
  left_join(study_n_df,by='study_id') %>% 
  rename(df=parameter)
write.csv(res_df,'adjustment_mode_ttest0.csv',row.names=F)

res_df=df %>% 
  select(orig_es,bayes_mean_es,sample_size,roi,study_id) %>%
  mutate(dist_es=abs(orig_es-bayes_mean_es)) %>% 
  group_by(study_id) %>% 
  nest() %>% 
  mutate(t_model=map(data,~t.test(.x$dist_es,mu=0,alternative='greater'))) %>%  #
  mutate(tidy_db =  map(t_model,tidy)) %>% # 
  unnest(tidy_db) %>% 
  select(study_id,estimate, statistic, p.value, parameter) %>% 
  left_join(study_n_df,by='study_id') %>% 
  rename(df=parameter)
write.csv(res_df,'adjustment_mean_ttest0.csv',row.names=F)


