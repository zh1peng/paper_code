#===============================================
#1. Define functions for cohen's d
d <- function(df, t, n1, n2) {
  d <- t*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
  return(d)}

se.d <-function(d,n1,n2){
  se<-sqrt((n1+n2)/(n1*n2)+(d^2)/(2*(n1+n2-2)))
  names(se)<-"se for d"
  return(se)}

ci.d <- function(d, SE){
  ci.d <- c((d-1.96*SE),(d+1.96*SE))
  names(ci.d) <- c("95% CI lower", "95% CI upper")
  return(ci.d)
}


#===============================================
#2 do group comp with site as random effects
do_group_comp <- function(df,if_meanCT=F){
  df$meanCT=rowMeans(df[grep("^L_|^R_", colnames(df))])
  base_fomular='region_data~Group+Sex + Age + ICV'
  if (if_meanCT){
    base_fomular=paste(base_fomular,'meanCT',sep='+')
  }
  
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(df)[grep('^L_|^R_',colnames(df))]){
    df$region_data=df[,region_i]
    lme_model=lme(as.formula(base_fomular),random=~1|Site,data=df,
                  method='REML',control = lmeControl(opt = "optim"))
    lme_model$call$fixed <- eval(lme_model$call$fixed) # add fomular information for the subsequent checks
    t_g[[region_i]]=summary(lme_model)$tTable["Group1",'t-value']
    p_g[[region_i]]=summary(lme_model)$tTable["Group1",'p-value']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lme_model)$tTable["Group1","DF"],
                         summary(lme_model)$tTable["Group1","t-value"],
                         as.numeric(table(df$Group)["0"]),
                         as.numeric(table(df$Group)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$Group)["0"]),
                             as.numeric(table(df$Group)["1"]))
    ## need to trace(performance:::check_normality.default, edit=TRUE)
    # to modify the residualize method used for the LME model as stats::rstandard doesn't work
    plot.n[[region_i]]=plot(check_normality(lme_model))+ggtitle(region_i)
    plot.v[[region_i]]=plot(check_heteroscedasticity(lme_model))+ggtitle(region_i)
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  rownames(results)=colnames(data.frame(t_g))
  
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}

#===============================================
#3 do group comp with SubjID as random effects (longitudinal analysis)
do_dev_comp <- function(df,if_meanCT=F, if_euler_score=T, performance=F){
  df$meanCT=rowMeans(df[grep("^L_|^R_", colnames(df))])
  base_fomular='region_data~time+Sex+mod_PDS+ICV+Site+time_diff'
  if (if_meanCT){
    base_fomular=paste(base_fomular,'meanCT',sep='+')
  }
  if (if_euler_score){
    base_fomular=paste(base_fomular,'qc_euler_score',sep='+')
  }
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(df)[grep('^L_|^R_',colnames(df))]){
    df$region_data=df[,region_i]
    lme_model=lme(as.formula(base_fomular),random=~1|SubjID,data=df,
                  method='REML',control = lmeControl(opt = "optim"))
    lme_model$call$fixed <- eval(lme_model$call$fixed) # add fomular information for the subsequent checks
    t_g[[region_i]]=summary(lme_model)$tTable["time1",'t-value']
    p_g[[region_i]]=summary(lme_model)$tTable["time1",'p-value']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lme_model)$tTable["time1","DF"],
                         summary(lme_model)$tTable["time1","t-value"],
                         as.numeric(table(df$time)["0"]),
                         as.numeric(table(df$time)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$time)["0"]),
                             as.numeric(table(df$time)["1"]))
    plot.n[[region_i]]=plot(check_normality(lme_model))+ggtitle(region_i)
    plot.v[[region_i]]=plot(check_heteroscedasticity(lme_model))+ggtitle(region_i)
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  rownames(results)=colnames(data.frame(t_g))
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}


#===============================================
#4 do group comp with resid data
do_group_comp_with_resid <- function(df,...){
  df$meanCT=rowMeans(df[grep("^L_|^R_", colnames(df))])
  resid_df=resid_CT(df,...)
  resid_df$Group=df$Group
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(resid_df)[grep('^resid_L_|^resid_R_',colnames(resid_df))]){
    resid_df$region_data=resid_df[,region_i]
    lm_model=lm(region_data~Group,data=resid_df)
    t_g[[region_i]]=summary(lm_model)$coefficients["Group1",'t value']
    p_g[[region_i]]=summary(lm_model)$coefficients["Group1",'Pr(>|t|)']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lm_model)$df[2],
                         summary(lm_model)$coefficients["Group1",'t value'],
                         as.numeric(table(df$Group)["0"]),
                         as.numeric(table(df$Group)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$Group)["0"]),
                             as.numeric(table(df$Group)["1"]))
    plot.n[[region_i]]=plot(check_normality(lm_model))+ggtitle(region_i)
    plot.v[[region_i]]= plot(check_homogeneity(lm_model))+ggtitle(region_i)
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}



#===============================================
#5 do group comp with pc1 removed data
do_group_comp_with_pc1_rm <- function(df,...){
  resid_df=resid_CT(df,...)
  pca_data=as.matrix(resid_df)
  pca=prcomp(pca_data, center = T, scale. = T) 
  pca=remove_pc1(pca)
  df2lm=pca$pc1_rm_data
  df2lm$Group=df$Group
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(df2lm)[grep('^pc1_rm_resid_L_|^pc1_rm_resid_R_',colnames(df2lm))]){
    df2lm$region_data=df2lm[,region_i]
    lm_model=lm(region_data~Group,data=df2lm)
    t_g[[region_i]]=summary(lm_model)$coefficients["Group1",'t value']
    p_g[[region_i]]=summary(lm_model)$coefficients["Group1",'Pr(>|t|)']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lm_model)$df[2],
                         summary(lm_model)$coefficients["Group1",'t value'],
                         as.numeric(table(df$Group)["0"]),
                         as.numeric(table(df$Group)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$Group)["0"]),
                             as.numeric(table(df$Group)["1"]))
    plot.n[[region_i]]=plot(check_normality(lm_model))+ggtitle(region_i)
    plot.v[[region_i]]=plot(check_homogeneity(lm_model))+ggtitle(region_i)+ggtitle(region_i)
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}

#===============================================
#6 do long comp with resid data
do_long_comp_with_resid <- function(df,f2resid){
  resid_df=resid_CT_by_formula(df,f2resid)
  resid_df$time=df$time
  resid_df$SubjID=df$SubjID
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(resid_df)[grep('^resid_L_|^resid_R_',colnames(resid_df))]){
    resid_df$region_data=resid_df[,region_i]
    lme_model=lme(region_data~time, random=~1|SubjID, data = resid_df,
                  method='REML',control = lmeControl(opt = "optim"))
    t_g[[region_i]]=summary(lme_model)$tTable["time1",'t-value']
    p_g[[region_i]]=summary(lme_model)$tTable["time1",'p-value']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lme_model)$tTable["time1","DF"],
                         summary(lme_model)$tTable["time1","t-value"],
                         as.numeric(table(df$time)["0"]),
                         as.numeric(table(df$time)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$time)["0"]),
                             as.numeric(table(df$time)["1"]))
    plot.n[[region_i]]=plot(check_normality(lme_model))+ggtitle(region_i)
    plot.v[[region_i]]=plot(check_heteroscedasticity(lme_model))+ggtitle(region_i)
    
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  rownames(results)=colnames(data.frame(t_g))
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}


#===============================================
#7 do long comp with pc1 removed
do_long_comp_with_pc1_rm <- function(df,f2resid){
  resid_df=resid_CT_by_formula(df,f2resid)
  pca_data=as.matrix(resid_df)
  pca=prcomp(pca_data, center = T, scale. = T) 
  pca=remove_pc1(pca)
  df2lme=pca$pc1_rm_data
  df2lme$time=df$time
  df2lme$SubjID=df$SubjID
  t_g=list()
  p_g=list()
  e_g=list()
  se_g=list()
  plot.n=list()
  plot.v=list()
  for (region_i in colnames(df2lme)[grep('^pc1_rm_resid_L_|^pc1_rm_resid_R_',colnames(df2lme))]){
    df2lme$region_data=df2lme[,region_i]
    lme_model=lme(region_data~time, random=~1|SubjID, data = df2lme,
                  method='REML',control = lmeControl(opt = "optim"))
    t_g[[region_i]]=summary(lme_model)$tTable["time1",'t-value']
    p_g[[region_i]]=summary(lme_model)$tTable["time1",'p-value']
    # e_g[[region_i]]=standardize_parameters(lme_model)$Std_Coefficient[2]
    e_g[[region_i]] <- d(summary(lme_model)$tTable["time1","DF"],
                         summary(lme_model)$tTable["time1","t-value"],
                         as.numeric(table(df$time)["0"]),
                         as.numeric(table(df$time)["1"]))
    se_g[[region_i]] <- se.d(e_g[[region_i]],
                             as.numeric(table(df$time)["0"]),
                             as.numeric(table(df$time)["1"]))
    plot.n[[region_i]]=plot(check_normality(lme_model))+ggtitle(region_i)
    plot.v[[region_i]]=plot(check_heteroscedasticity(lme_model))+ggtitle(region_i)
    
  }
  
  results=data.frame(t_vals=unlist(t_g),
                     p_vals=unlist(p_g),
                     d_vals=unlist(e_g),
                     se_vals=unlist(se_g)) %>% 
    mutate(ci95_low=d_vals-1.96*se_vals,
           ci95_high=d_vals+1.96*se_vals,
           fdr_p=p.adjust(p_vals, method="fdr"),
           if_sig=fdr_p<0.05)
  rownames(results)=colnames(data.frame(t_g))
  res=list(results=results,
           plot_norm=plot.n,
           plot_variance=plot.v)
  return(res)
}



#===============================================
#8 plot results
# if_sig: if plot effect size maps with p<0.05
plot_group_comp <- function(results,limits=c(-0.2,0.2),showlegend=F, if_sig=F){
  if (if_sig){
    results= results%>%filter(results$fdr_p<0.05) 
  }
  p=results %>% 
    mutate(label=rownames(results),
           label=sub('.*L_','lh_',label),
           label=sub('.*R_','rh_',label),
           label=sub('_thickavg','',label))%>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = d_vals)) +
    geom_sf(show.legend = showlegend)+
    scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=limits,aes(title='ES'))+
    theme_void()
  
  p
  return(p)
}







