library(dplyr)
cbind.fill <- function(...){ # can use list instead
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


setwd('F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model')
filenames=c('L_ct', 'L_sa', 'L_sv', 'R_ct', 'R_sa', 'R_sv')

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



# 1. Export table
key_df=read.csv('regions.csv')
orig_label=as.character(key_df$orignal_label)
full_label=as.character(key_df$full_label)
ggseg_label=as.character(key_df$ggseg_label)

hem_mu=all_ROImean[,colnames(all_ROImode)[grep('L_',colnames(all_ROImode))]]
hem_ess=all_ROIess[,colnames(all_ROIess)[grep('L_',colnames(all_ROIess))]]
hem_lower=all_ROIlower[,colnames(all_ROIlower)[grep('L_',colnames(all_ROIlower))]]
hem_upper=all_ROIupper[,colnames(all_ROIupper)[grep('L_',colnames(all_ROIupper))]]
hem_rhat=all_ROIrhat[,colnames(all_ROIrhat)[grep('L_',colnames(all_ROIrhat))]]

hem_df=data.frame(matrix(nrow=length(colnames(hem_mu)[grep("L_", colnames(hem_mu))]), ncol=0),stringsAsFactors = F)
hem_df[,"region"]= sub("L_", "", colnames(hem_mu)[grep("L_", colnames(hem_mu))])

hem_df[,'mu']=round(as.numeric(hem_mu['mu',]),digits=3)

hem_df[,'lower']=round(as.numeric(hem_lower['mu',]),digits=3)
hem_df[,'upper']=round(as.numeric(hem_upper['mu',]),digits=3)

hem_df[,'ess']=round(as.numeric(hem_ess['mu',]),digits=3)
hem_df[,'rhat']=round(as.numeric(hem_rhat['mu',]),digits=3)

hem_df[,'region_full']=stringr::str_replace_all(hem_df[,"region"],
                                                setNames(full_label,paste0(orig_label,'$')))
write.csv(hem_df,'F:/Google Drive/post-doc/Bayesian_Project/new_model/brain_and_ridge_plots/left_M_table.csv')


hem_mu=all_ROImean[,colnames(all_ROImode)[grep('R_',colnames(all_ROImode))]]
hem_ess=all_ROIess[,colnames(all_ROIess)[grep('R_',colnames(all_ROIess))]]
hem_lower=all_ROIlower[,colnames(all_ROIlower)[grep('R_',colnames(all_ROIlower))]]
hem_upper=all_ROIupper[,colnames(all_ROIupper)[grep('R_',colnames(all_ROIupper))]]
hem_rhat=all_ROIrhat[,colnames(all_ROIrhat)[grep('R_',colnames(all_ROIrhat))]]
hem_df=data.frame(matrix(nrow=length(colnames(hem_mu)[grep("R_", colnames(hem_mu))]), ncol=0),stringsAsFactors = F)
hem_df[,"region"]= sub("R_", "", colnames(hem_mu)[grep("R_", colnames(hem_mu))])
hem_df[,'mu']=round(as.numeric(hem_mu['mu',]),digits=3)
hem_df[,'lower']=round(as.numeric(hem_lower['mu',]),digits=3)
hem_df[,'upper']=round(as.numeric(hem_upper['mu',]),digits=3)
hem_df[,'ess']=round(as.numeric(hem_ess['mu',]),digits=3)
hem_df[,'rhat']=round(as.numeric(hem_rhat['mu',]),digits=3)
hem_df[,'region_full']=stringr::str_replace_all(hem_df[,"region"],
                                              setNames(full_label,paste0(orig_label,'$')))
write.csv(hem_df,'F:/Google Drive/post-doc/Bayesian_Project/new_model/brain_and_ridge_plots/right_M_table.csv')








