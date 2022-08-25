# Code for the analysis 
* Analysis 1. PCA across datasets
* Analysis 2. Case-control comparisons
* Analysis 3. Neurodevelopmental effects
* Analysis 4. Gene expression analysis

## Analysis 1-3 used R 3.6.2 
SessionInfo:
```
R version 3.6.2 (2019-12-12)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] see_0.7.2          performance_0.9.2  matrixStats_0.58.0 Hmisc_4.3-1        Formula_1.2-4      survival_3.1-8    
 [7] lattice_0.20-38    corrplot_0.84      nlme_3.1-143       gridExtra_2.3      ggseg_1.6.02       ggplot2_3.3.6     
[13] tidyr_1.1.3        dplyr_1.0.5       

loaded via a namespace (and not attached):
 [1] xfun_0.22           tidyselect_1.1.0    purrr_0.3.4         sf_0.9-8            splines_3.6.2      
 [6] colorspace_2.0-0    vctrs_0.3.7         generics_0.1.0      htmltools_0.5.1.1   base64enc_0.1-3    
[11] utf8_1.2.1          rlang_0.4.12        e1071_1.7-6         pillar_1.6.0        foreign_0.8-74     
[16] glue_1.4.2          withr_2.4.2         DBI_1.1.1           RColorBrewer_1.1-2  jpeg_0.1-8.1       
[21] lifecycle_1.0.0     plyr_1.8.6          stringr_1.4.0       munsell_0.5.0       gtable_0.3.0       
[26] htmlwidgets_1.5.1   knitr_1.33          latticeExtra_0.6-29 class_7.3-17        fansi_0.4.2        
[31] htmlTable_1.13.3    Rcpp_1.0.6          acepack_1.4.1       KernSmooth_2.23-16  backports_1.2.1    
[36] scales_1.1.1        classInt_0.4-3      checkmate_2.0.0     digest_0.6.27       png_0.1-7          
[41] stringi_1.5.3       insight_0.18.2      tools_3.6.2         magrittr_2.0.1      proxy_0.4-25       
[46] tibble_3.1.1        cluster_2.1.0       crayon_1.4.1        pkgconfig_2.0.3     ellipsis_0.3.1     
[51] Matrix_1.2-18       data.table_1.13.2   rstudioapi_0.13     assertthat_0.2.1    R6_2.5.0           
[56] rpart_4.1-15        units_0.7-1         nnet_7.3-12         compiler_3.6.2 
```


```R

# load things
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg) # brain plot
library(gridExtra)
library(grid)
library(nlme)
library(corrplot) # make corrplot
library(Hmisc) # rcorr
library(matrixStats) # spin test
library(performance)
# need to trace(performance:::check_normality.default, edit=TRUE)
# to modify the residualize method used for the LME model as stats::rstandard doesn't work for lme
# if (class(x)=='lm'){
#   p.val <- .check_normality(stats::rstandard(x), x)
# }else if (class(x)=='lme'){
#   p.val=.check_normality(residuals(x, type = "normalized"),x)
# }
library(see)

data_path='F:/Google Drive/post-doc/p_factor'
code_path='F:/Google Drive/post-doc/p_factor/clean_code'
source(paste(code_path,'pca_functions.R',sep = '/'))
source(paste(code_path,'group_comp_functions.R',sep = '/'))
source(paste(code_path,'correlation_functions.R',sep = '/'))
source(paste(code_path,'external_functions.R',sep = '/'))

# create and save perm.id 
# it will be used for all spin tests
# perm.id=get_68perm.id(10000)
# saveRDS(perm.id,paste(code_path,'68perm_id.rds',sep='/'))
perm.id=readRDS(paste(code_path,'68perm_id.rds',sep='/')) 


# Analysis 1: PCA analysis
result_path=file.path(data_path,'all_results/analysis1')
dir.create(result_path)

#1.1. do PCA==================================
pca.abcd_t1=do_pca(data_path, data_name = 'ABCD',if_pds = T,if_resid=T,if_meanCT = F)
pca.abcd_t2=do_pca(data_path, data_name='ABCD_FU',if_pds=T, if_resid=T, if_meanCT=F)
pca.imagen_t1=do_pca(data_path, data_name = 'IMAGEN_BSL',if_pds = T,if_resid=T,if_meanCT = F)
pca.imagen_t2=do_pca(data_path, data_name = 'IMAGEN_FU2',if_resid=T,if_meanCT = F)
pca.enigma_hc=do_pca(data_path, data_name = 'ENIGMA_HC',if_resid=T,if_meanCT = F)
pca.ukb=do_pca(data_path, data_name='UKB',if_resid=T,if_meanCT = F)
pca.merged=do_pca(data_path, data_name = 'MERGED',if_resid=T,if_meanCT = F)

# Test if the order of the variables matters [doesn't matter]
# pca.merged.test=do_pca(data_path, data_name = 'MERGEDRANDOM',if_resid=T,if_meanCT = F)
# sort(pca.merged.test$stand_loadings[,'PC1'])
# sort(pca.merged$stand_loadings[,'PC1'])

#1.2. plot all PC1==========================
pca_list=list(ABCD_T1=pca.abcd_t1,
              ABCD_T2=pca.abcd_t2,
              IMAGEN_T1=pca.imagen_t1,
              IMAGEN_T2=pca.imagen_t2,
              UKB=pca.ukb,
              ENIGMA_HC=pca.enigma_hc,
              Merged=pca.merged)

# make plot for each of them
pc2plot=1
for (data_i in names(pca_list)){
  pca.obj=eval(parse(text=sprintf('pca_list$%s',data_i))) # use this way to get the prcomp data format 
  p=plot_pca_brain(pca.obj,pc2plot = pc2plot)+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  ggsave(sprintf('%s/%s_pc%s.png',result_path, data_i[1],pc2plot),p,width = 6, height = 4)
}

# test a different color scheme
for (data_i in names(pca_list)){
  pca.obj=eval(parse(text=sprintf('pca_list$%s',data_i))) # use this way to get the prcomp data format 
  p=plot_pca_brain2(pca.obj,pc2plot = pc2plot)+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  ggsave(sprintf('%s/%s_pc%s_v1.png',result_path, data_i[1],pc2plot),p,width = 6, height = 4)
}
# get a colorbar from one plot
p_colorbar=plot_pca_brain2(pca.obj,pc2plot = pc2plot,T)
ggsave(sprintf('%s/%s_pc%s_v1_withcolorbar.png',result_path, data_i[1],pc2plot),p_colorbar,width = 6, height = 4)


#make pca_plot
for (data_i in names(pca_list)){
  pca.obj=eval(parse(text=sprintf('pca_list$%s',data_i))) # use this way to get the prcomp data format 
  p=plot_pca(pca.obj,title2add = data_i)
  ggsave(sprintf('%s/%s_pcaplot_pc.png',result_path, data_i[1]),p,width = 6, height = 10)
}



#1.3. do correlation for PCs
pc2plot='PC1'
df2cor=cbind(pca.ukb$stand_loadings[,pc2plot], 
             pca.enigma_hc$stand_loadings[,pc2plot],
             pca.imagen_t2$stand_loadings[,pc2plot],
             pca.abcd_t1$stand_loadings[,pc2plot],  
             pca.abcd_t2$stand_loadings[,pc2plot],
             pca.imagen_t1$stand_loadings[,pc2plot])
colnames(df2cor) <- c('UKB',
                      'ENIGMA-CTRL',
                      'IMAGEN-T2',
                      'ABCD-T1',
                      'ABCD-T2',
                      'IMAGEN-T1')
cor.list=rcorr(as.matrix(df2cor),type='pearson')
cor.list$p.spin=caculate_spin_p(df2cor,perm.id)

# save correlation plot
png(file=sprintf("%s/%s_cor_mat_spin_sig.png",result_path,pc2plot),
    width=800, height=800,res = 150)
cor_plot(cor.list,type='upper')
dev.off()

# save correlation results
t2save=cor_t2save(cor.list)
write.csv(t2save,sprintf('%s/%s_cor_mat_spin.csv',result_path,pc2plot))

#1.4. report loadings of the PCs; VAF; eigvalue
df2report=list()
for (data_i in names(pca_list)){
  pca.obj=eval(parse(text=sprintf('pca_list$%s',data_i)))
  sloadings=pca.obj$stand_loadings
  eigvalue=t(as.data.frame(pca.obj$sdev^2))
  colnames(eigvalue)=colnames(sloadings)
  rownames(eigvalue)='Eigenvalue'
  VAF=t(as.data.frame(pca.obj$sdev^2/sum(pca.obj$sdev^2)))*100
  colnames(VAF)=colnames(sloadings)
  rownames(VAF)='Variance accounted for (VAF %)'
  df2report[[data_i]]=rbind(sloadings,eigvalue,VAF)
}
df2save=do.call('cbind',df2report)
write.csv(df2report,sprintf('%s/all_PCA_results.csv',result_path))


#1.5. plot VAF for each dataset
# This is not a generic function, the dataset names were hardcoded in the function
p=plot_pca_list_VAF(pca_list,pc2plot=8,type = 'bar')
ggsave(sprintf('%s/VAF_plot.png',result_path),p,width = 8, height = 6)


# Analysis 2: Case-control Comparisons
result_path=file.path(data_path,'all_results/analysis2')
dir.create(result_path)

# 2.1. group comparison in ENIGMA
data_name='ENIGMA_CASE_CONTROL'
g_comp.enigma=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.enigma=preprocess_df(g_comp.enigma,Group = T)
result.enigma=do_group_comp(g_comp.enigma,if_meanCT = F)
write.csv(result.enigma$results,sprintf('%s/enigma_diff_results.csv',result_path))
# visual check normality
marrangeGrob(grobs=result.enigma$plot_norm, nrow=4, ncol=4)
# visual check heteroscedasticity
marrangeGrob(grobs=result.enigma$plot_variance, nrow=4, ncol=4)

# 2.2. group comparison in UKB Group (defined by AUDIT score)
data_name='UKB_L1_H19'
g_comp.ukb=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.ukb=preprocess_df(g_comp.ukb,Group = T)
write.csv(result.ukb$results,sprintf('%s/ukb_diff_results.csv',result_path))
# visual check normality
marrangeGrob(grobs=result.ukb$plot_norm, nrow=4, ncol=4)
# visual check heteroscedasticity
marrangeGrob(grobs=result.ukb$plot_variance, nrow=4, ncol=4)


# 2.3. plot case-control comparison effect size
# Note: effect size from group comparison was copied to this file
es_df=read.csv(sprintf('%s/all_case_control_es.csv',code_path)) # read case-control effect size
for (idx in c(2:length(colnames(es_df)))){
  col2plot=colnames(es_df)[idx]
  p=es_df%>%
    select(label,col2plot)%>%
    mutate(
      label=sub('.*L_','lh_',label),
      label=sub('.*R_','rh_',label),
      label=sub('_thickavg','',label))%>%
    rename('ES'=col2plot)%>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = ES)) +
    geom_sf(show.legend = F) +
    scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=c(-0.6,0.6),aes(title="ES"))+
    theme_void()
  ggsave(sprintf('%s/ES_%s.png',result_path, col2plot),p,width = 6, height = 4)
}



# 2.4. do ES correlation    
df2cor=cbind(es_df$ENIGMA_DIFF,
             es_df$UKB_DIFF,
             es_df$BD,
             es_df$MDD,
             es_df$OCD,
             es_df$SCZ,
             es_df$EPI,
             es_df$CHR,
             pca.merged$stand_loadings[,'PC1'])

colnames(df2cor) <- c('ALCgc (ENIGMA)',
                      'ALCgc (UKB)',
                      'BD',
                      'MDD',
                      'OCD',
                      'SCZ',
                      'EPI',
                      'CHR',
                      'Combined-PC1')

cor.list=rcorr(as.matrix(df2cor),type='pearson')
cor.list$p.spin=caculate_spin_p(df2cor,perm.id)
png(file=sprintf("%s/es_cor_mat_spin_sig.png", result_path),
    width=1200, height=1200,res = 150)
cor_plot(cor.list,type='upper')
dev.off()
t2save=cor_t2save(cor.list)
write.csv(t2save,sprintf('%s/es_cor_mat_spin.csv',result_path))


# 2.5. compare results with and without removing PC1

## group comp with resid
data_name='ENIGMA_CASE_CONTROL'
g_comp.enigma=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.enigma=preprocess_df(g_comp.enigma,Group = T)
result_resid.enigma=do_group_comp_with_resid(g_comp.enigma,if_meanCT = F)
# visual check on normality
marrangeGrob(grobs=result_resid.enigma$plot_norm, nrow=4, ncol=4)
# visual check on homogeneity
marrangeGrob(grobs=result_resid.enigma$plot_variance, nrow=4, ncol=4)


## group comp with pc1 removed
data_name='ENIGMA_CASE_CONTROL'
g_comp.enigma=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.enigma=preprocess_df(g_comp.enigma,Group = T)
result_resid_pc1rm.enigma=do_group_comp_with_pc1_rm(g_comp.enigma,if_meanCT = F)
# visual check on normality
marrangeGrob(grobs=result_resid_pc1rm.enigma$plot_norm, nrow=4, ncol=4)
# visual check on homogeneity
marrangeGrob(grobs=result_resid_pc1rm.enigma$plot_variance, nrow=4, ncol=4)

# DO same thing for UKB data
## group comp with resid
data_name='UKB_L1_H19'
g_comp.ukb=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.ukb=preprocess_df(g_comp.ukb,Group = T)
result_resid.ukb=do_group_comp_with_resid(g_comp.ukb,if_meanCT = F)
# visual check on normality
marrangeGrob(grobs=result_resid.ukb$plot_norm, nrow=4, ncol=4)
# visual check on homogeneity
marrangeGrob(grobs=result_resid.ukb$plot_variance, nrow=4, ncol=4)


## group comp with pc1 removed
data_name='UKB_L1_H19'
g_comp.ukb=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name))
g_comp.ukb=preprocess_df(g_comp.ukb,Group = T)
result_resid_pc1rm.ukb=do_group_comp_with_pc1_rm(g_comp.ukb,if_meanCT = F)
# visual check on normality
marrangeGrob(grobs=result_resid.ukb$plot_norm, nrow=4, ncol=4)
# visual check on homogeneity
marrangeGrob(grobs=result_resid.ukb$plot_variance, nrow=4, ncol=4)

result_list=list(ENIGMA_normal=result_resid.enigma$results,
                 ENIGMA_rm_PC1=result_resid_pc1rm.enigma$results,
                 UKB_normal=result_resid.ukb$results,
                 UKB_rm_PC1=result_resid_pc1rm.ukb$results)

# save plot and results
for (result_i in names(result_list)){
  df2use=eval(parse(text=sprintf('result_list$%s',result_i)))
  p1=plot_group_comp(df2use, limits = c(-0.6, 0.6), if_sig = F, showlegend = F)
  p2=plot_group_comp(df2use, limits = c(-0.6, 0.6), if_sig = T, showlegend = F)
  write.csv(df2use, sprintf('%s/supp_analysis_%s.csv',result_path,result_i))
  ggsave(sprintf('%s/supp_analysis_%s_sig.png',result_path,result_i),p2,width = 6,height = 4)
  ggsave(sprintf('%s/supp_analysis_%s_all.png',result_path,result_i),p1,width = 6,height = 4)
}


# 2.6. Compare es of ukb/enigma to Combined-PC1 without control participants from ENIGMA
pca.merged.noHC=do_pca(data_path, data_name = 'MERGED_noHC',if_resid=T,if_meanCT = F)
df2cor=cbind(es_df$ENIGMA_DIFF,
             es_df$UKB_DIFF,
             es_df$BD,
             es_df$MDD,
             es_df$OCD,
             es_df$SCZ,
             es_df$EPI,
             es_df$CHR,
             pca.merged.noHC$stand_loadings[,'PC1'])

colnames(df2cor) <- c('ALCgc(ENIGMA)',
                      'ALCgc(UKB)',
                      'BD',
                      'MDD',
                      'OCD',
                      'SCZ',
                      'EPI',
                      'CHR',
                      'Combined-PC1(without control participants)')

cor.list=rcorr(as.matrix(df2cor),type='pearson')
cor.list$p.spin=caculate_spin_p(df2cor,perm.id)
png(file=sprintf("%s/supp_es_cor_mat_spin_sig.png", result_path),
    width=1200, height=1200,res = 150)
cor_plot(cor.list,type='upper')
dev.off()
t2save=cor_t2save(cor.list)
write.csv(t2save,sprintf('%s/supp_es_cor_mat_spin.csv',result_path))


# Analysis 3: Neurodevelopmental effects
result_path=file.path(data_path,'all_results/analysis3')
dir.create(result_path)

#3.1. compare BSL to FU2
## 3.1.1. read data and add euler number
data_name='IMAGEN_BSL'
df.imagen.t1=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name)) 
df.imagen.t1=preprocess_df(df.imagen.t1)%>% mutate(time=0) 
qc_score=read.csv(sprintf('%s/%s/euler.tsv',data_path,data_name),sep = '\t') %>% 
  rename(SubjID=orig.nofix) %>% 
  rename (qc_euler_score=mean) %>% 
  select(SubjID, qc_euler_score)

df.imagen.t1=plyr::join(df.imagen.t1,qc_score,'SubjID')%>%
  filter(complete.cases(.))
df.imagen.t1$qc_euler_score=scale(df.imagen.t1$qc_euler_score)

data_name='IMAGEN_FU2'
df.imagen.t2=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name)) 
df.imagen.t2=preprocess_df(df.imagen.t2)%>% mutate(time=1)
qc_score=read.csv(sprintf('%s/%s/euler.tsv',data_path,data_name),sep = '\t') %>% 
  rename(SubjID=orig.nofix) %>% 
  rename (qc_euler_score=mean) %>% 
  select(SubjID, qc_euler_score)
df.imagen.t2=plyr::join(df.imagen.t2,qc_score,'SubjID')%>%
  filter(complete.cases(.))
df.imagen.t2$qc_euler_score=scale(df.imagen.t2$qc_euler_score)


## add sampling deviation from mean age at T1 and T2 for each participant
df.imagen.t1$time_diff=df.imagen.t1$Age-mean(df.imagen.t1$Age)
df.imagen.t2$time_diff=df.imagen.t2$Age-mean(df.imagen.t2$Age)
df.imagen.all=rbind(df.imagen.t1,df.imagen.t2)
df.imagen.all$SubjID=as.factor(df.imagen.all$SubjID)
df.imagen.all$time=as.factor(df.imagen.all$time)
df.imagen.all$qc_euler_score=as.numeric(df.imagen.all$qc_euler_score)
df.imagen.all$time_diff=as.numeric(df.imagen.all$time_diff)
result.imagen.dev=do_dev_comp(df.imagen.all,if_euler_score=T)
write.csv(result.imagen.dev$results,sprintf('%s/IMAGEN_dev_diff_results.csv',result_path))

# visual check normality
marrangeGrob(grobs=result.imagen.dev$plot_norm, nrow=4, ncol=4)
# visual check heteroscedasticity
marrangeGrob(grobs=result.imagen.dev$plot_variance, nrow=4, ncol=4)

# 3.2. plot effect size maps [this has to be separate as it is in a very different scale]
p=plot_group_comp(result.imagen.dev$results, limits=c(-2,2), showlegend = F)
ggsave(sprintf('%s/IMAGEN_diff.png',result_path),p,width = 6, height = 4)



# 3.3. ABCD T1 vs.T2
data_name='ABCD'
df.abcd.t1=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name)) 
df.abcd.t1=preprocess_df(df.abcd.t1)%>% mutate(time=0) %>% select(-BSL_pds)
data_name='ABCD_FU'
df.abcd.t2=read.csv(sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name)) 
df.abcd.t2=preprocess_df(df.abcd.t2)%>% mutate(time=1) %>% select(-FU_pds)

df.abcd.t1$time_diff=df.abcd.t1$Age-mean(df.abcd.t1$Age)
df.abcd.t2$time_diff=df.abcd.t2$Age-mean(df.abcd.t2$Age)
df.abcd.all=rbind(df.abcd.t1,df.abcd.t2)
df.abcd.all$SubjID=as.factor(df.abcd.all$SubjID)
df.abcd.all$time=as.factor(df.abcd.all$time)
df.abcd.all$time_diff=as.numeric(df.abcd.all$time_diff)
result.abcd.dev=do_dev_comp(df.abcd.all,if_euler_score=F)
write.csv(result.abcd.dev$results,sprintf('%s/ABCD_dev_diff_results.csv',result_path))
# visual check normality
marrangeGrob(grobs=result.abcd.dev$plot_norm, nrow=4, ncol=4)
# visual check heteroscedasticity
marrangeGrob(grobs=result.abcd.dev$plot_variance, nrow=4, ncol=4)


# 3.4. plot effect size maps
p=plot_group_comp(result.abcd.dev$results, limits=c(-2,2), showlegend = F)
ggsave(sprintf('%s/ABCD_diff.png',result_path),p,width = 6, height = 4)


# 3.5. plot results from Lifespan paper
es_df=read.csv(sprintf('%s/all_dev_es.csv',code_path))
for (idx in c(2:length(colnames(es_df)))){
  col2plot=colnames(es_df)[idx]
  p=es_df%>%
    select(label,col2plot)%>%
    mutate(
      label=sub('.*L_','lh_',label),
      label=sub('.*R_','rh_',label),
      label=sub('_thickavg','',label))%>%
    rename('ES'=col2plot)%>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = ES)) +
    geom_sf(show.legend = F) +
    scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=c(-0.6,0.6),aes(title="ES"))+
    theme_void()
  ggsave(sprintf('%s/ES_%s.png',result_path, col2plot),p,width = 6, height = 4)
}


# 3.6. do correlation analysis
df2cor=cbind(result.abcd.dev$results$d_vals,
             result.imagen.dev$results$d_vals,
             es_df$young_age_r,
             es_df$middle_age_r,
             es_df$old_age_r,
             es_df$LIFE_AGE,
             pca.merged$stand_loadings[,'PC1'])
colnames(df2cor) <- c('ABCD 10->12 years',
                      'IMAGEN 14->19 years',
                      'Correlation (3-29 years)',
                      'Correlation (30-59 years)',
                      'Correlation (60-90 years)',
                      'Variance (3-90 years)',
                      'Combined-PC1')
cor.list=rcorr(as.matrix(df2cor),type='pearson')
cor.list$p.spin=caculate_spin_p(df2cor,perm.id)
png(file=sprintf('%s/dev_cor_mat_spin_sig.png',result_path),
    width=1200, height=1200,res = 150)
cor_plot(cor.list,type='upper')
dev.off()
t2save=cor_t2save(cor.list)
write.csv(t2save,sprintf('%s/dev_cor_mat_spin.csv',result_path))


# 3.7 IMAGEN.dev vs. combined.PC1 without IMAGEN-T2 participants
pca.merged.noT2=do_pca(data_path, data_name = 'MERGED_noT2',if_resid=T,if_meanCT = F)
df2cor=cbind(result.abcd.dev$results$d_vals,
             result.imagen.dev$results$d_vals,
             es_df$young_age_r,
             es_df$middle_age_r,
             es_df$old_age_r,
             es_df$LIFE_AGE,
             pca.merged.noT2$stand_loadings[,'PC1'])
colnames(df2cor) <- c('ABCD 10->12 years',
                      'IMAGEN 14->19 years',
                      'Correlation (3-29 years)',
                      'Correlation (30-59 years)',
                      'Correlation (60-90 years)',
                      'Variance (3-90 years)',
                      'Combined-PC1(without IMAGEN-T2)')
cor.list=rcorr(as.matrix(df2cor),type='pearson')
cor.list$p.spin=caculate_spin_p(df2cor,perm.id)
png(file=sprintf('%s/supp_analysis_dev_cor_mat_spin_sig.png',result_path),
    width=1200, height=1200,res = 150)
cor_plot(cor.list,type='upper')
dev.off()
t2save=cor_t2save(cor.list)
write.csv(t2save,sprintf('%s/supp_analysis_dev_cor_mat_spin.csv',result_path))

```

## Analysis 4 used R 4.1.0
SessionInfo:
```
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggplot2_3.3.5         GOSemSim_2.18.0       enrichplot_1.12.2     DOSE_3.18.1           org.Hs.eg.db_3.13.0  
 [6] AnnotationDbi_1.54.1  IRanges_2.26.0        S4Vectors_0.30.0      Biobase_2.52.0        BiocGenerics_0.38.0  
[11] clusterProfiler_4.0.2

loaded via a namespace (and not attached):
 [1] nlme_3.1-152           bitops_1.0-7           ggtree_3.0.2           bit64_4.0.5            RColorBrewer_1.1-3    
 [6] httr_1.4.2             GenomeInfoDb_1.28.1    tools_4.1.0            utf8_1.2.1             R6_2.5.0              
[11] lazyeval_0.2.2         DBI_1.1.1              colorspace_2.0-2       withr_2.5.0            tidyselect_1.1.1      
[16] gridExtra_2.3          bit_4.0.4              compiler_4.1.0         cli_3.0.1              scatterpie_0.1.6      
[21] shadowtext_0.0.8       scales_1.1.1           stringr_1.4.0          digest_0.6.27          XVector_0.32.0        
[26] pkgconfig_2.0.3        fastmap_1.1.0          rlang_1.0.4            rstudioapi_0.13        RSQLite_2.2.7         
[31] farver_2.1.0           generics_0.1.0         jsonlite_1.7.2         BiocParallel_1.26.1    dplyr_1.0.7           
[36] RCurl_1.98-1.3         magrittr_2.0.1         GO.db_3.13.0           GenomeInfoDbData_1.2.6 patchwork_1.1.1       
[41] Matrix_1.3-4           Rcpp_1.0.7             munsell_0.5.0          fansi_0.5.0            ape_5.5               
[46] viridis_0.6.1          lifecycle_1.0.1        stringi_1.6.2          ggraph_2.0.5           MASS_7.3-54           
[51] zlibbioc_1.38.0        plyr_1.8.6             qvalue_2.24.0          grid_4.1.0             blob_1.2.1            
[56] ggrepel_0.9.1          DO.db_2.9              crayon_1.4.1           lattice_0.20-44        graphlayouts_0.7.1    
[61] Biostrings_2.60.1      cowplot_1.1.1          splines_4.1.0          KEGGREST_1.32.0        pillar_1.6.1          
[66] fgsea_1.18.0           igraph_1.2.6           reshape2_1.4.4         fastmatch_1.1-0        glue_1.4.2            
[71] downloader_0.4         BiocManager_1.30.16    data.table_1.14.0      treeio_1.16.1          png_0.1-7             
[76] vctrs_0.3.8            tweenr_1.0.2           gtable_0.3.0           purrr_0.3.4            polyclip_1.10-0       
[81] tidyr_1.1.3            cachem_1.0.5           ggforce_0.3.3          tidygraph_1.2.0        tidytree_0.3.4        
[86] viridisLite_0.4.0      tibble_3.1.2           rvcheck_0.1.8          aplot_0.0.6            memoise_2.0.0         
[91] ellipsis_0.3.2        
```
```R
# 4. gene expression analysis
# need to change to R 4+ as clusterProfiler required
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(GOSemSim)
library(ggplot2)

result_path=file.path(data_path,'all_results/analysis4-ORA')
dir.create(result_path)

# 4.1. prepare data: merge pc1/gene/centroid
r2test=0.4 # 
gene_data=read.csv(sprintf('%s/allgenes_stable_r%s.csv',code_path,r2test))%>%rename_with(~paste0('gene_',.),-c('label'))
# PC1 pattern
pc1=pca.merged$stand_loadings %>% 
  mutate(label=rownames(.),
         label=sub('resid_L_','L_',label),
         label=sub('resid_R_','R_',label),
         label=sub('_thickavg','',label)) %>% 
  select(label,PC1) %>%
  rename(brain_PC1=PC1)

# coordinate data
lh_centroid_data=read.csv(sprintf('%s/lh_centroid.csv',code_path))
rh_centroid_data=read.csv(sprintf('%s/rh_centroid.csv',code_path))
lr_centroid=rbind(lh_centroid_data,rh_centroid_data)
colnames(lr_centroid)=c('label','centroid1','centroid2','centroid3')

# merge all data by label and remove regions with NA values
df_all=merge(merge(gene_data,pc1,by='label'),lr_centroid,by='label')%>%filter(complete.cases(.))


# 4.2. Find significant genes
gene_df=df_all %>% select(starts_with('gene_'))
# get centroid.l/centroid.r for 66 regions
centroid.l=df_all %>% filter(grepl('^L_',label)) %>%  select(starts_with('centroid'))
centroid.r=df_all %>% filter(grepl('^R_',label)) %>%  select(starts_with('centroid'))
brain_df=df_all %>% select(starts_with('brain'))

#perm.id=rotate.parcellation(as.matrix(centroid.l),as.matrix(centroid.r),nrot=10000)
#saveRDS(perm.id,paste(code_path,'66perm_id.rds',sep='/'))
perm.id=readRDS(paste(code_path,'66perm_id.rds',sep='/')) 

null_brain_df=build_null_brain_df(perm.id,brain_df)

true_r_vals=cor(gene_df,brain_df)
hist(true_r_vals,breaks = seq(min(true_r_vals), max(true_r_vals), length.out = 100),main = ' Distribution of true r')

null_r_vals=cor(gene_df,null_brain_df,method='pearson')
pos_null_dist=apply(null_r_vals,2,'max') # for each iteration find the max number
neg_null_dist=apply(null_r_vals,2,'min') # for each iteration find the min number

pos_cutoff=quantile(pos_null_dist,.95)
neg_cutoff=quantile(neg_null_dist,.05)

# calculate p value following the FWE procedure. 
fwe_spin_p=as.numeric()
for (n in c(1:length(true_r_vals))){
  fwe_spin_p[n]=caculate_p_from_null(pos_null_dist,neg_null_dist, true_r_vals[n])
}

# calculate a normal p spin value for each correlation
# p.spin <-matrix(nrow = ncol(gene_df), 1)
# for (i in c(1:ncol(gene_df))) {
#   print(i)
#   p.spin[i,1] <- perm.sphere.p(brain_df$brain_PC1,gene_df[,i],perm.id,corr.type='pearson')}

# calculate a normal p value without spin test
true_p_vals=outer(gene_df, brain_df, Vectorize(function(a, b) cor.test(a, b, method='pearson')$p.value))

all_results=data.frame(gene_name=sub('gene_','',colnames(gene_df)),
                       r_vals=true_r_vals[,1],
                       p_vals=true_p_vals[,1],
                       fwe_spin_p=fwe_spin_p) 

all_results$if_sig[all_results$r_vals>pos_cutoff & all_results$r_vals>0] = "Pos"
all_results$if_sig[all_results$r_vals<neg_cutoff & all_results$r_vals<0 ] = "Neg"
all_results$if_sig[all_results$r_vals>=neg_cutoff & all_results$r_vals <= pos_cutoff ] = "NS"
df2save=all_results %>% filter(if_sig!='NS') %>% select(-if_sig) %>% arrange(desc(r_vals))

write.csv(df2save, sprintf('%s/sig_gene_info_r%s.csv',result_path,r2test),row.names = F)

# 4.3. ORA Analysis
enricher_internal <- getFromNamespace('enricher_internal','clusterProfiler')

r2test=0.4
sig_genes_df=read.csv(sprintf('%s/sig_gene_info_r%s.csv',result_path,r2test))
all_genes_df=read.csv(sprintf('%s/allgenes_stable_r%s.csv',code_path,r2test))
sig_genes=sig_genes_df$gene_name

sig_gene_symbol=as.character(sig_genes)
bg_genes=colnames(all_genes_df)[-1]
bg_gene_symbol=as.character(bg_genes)
  
 # GO enrichment 
for (ont2use in c('BP','MF','CC')){
enrich.go=enrichGO(gene=sig_gene_symbol,
                  universe = bg_gene_symbol,
                  OrgDb=org.Hs.eg.db,
                  keyType="SYMBOL",
                  ont=ont2use,
                  pvalueCutoff = 1,
                  minGSSize = 10,
                  maxGSSize = 200)
ora_df=data.frame(enrich.go)
write.csv(ora_df,sprintf('%s/GO_%s_ORA_results_%s.csv',result_path,ont2use,r2test))
}

  
# ORA with development
for (pSI.threshold in c(0.01)){# 0.01 is used given the number of genes after intersection with gene expression data
  GS_file=sprintf('%s/pSI/pSI_only_cortex_table_%s.csv',code_path,pSI.threshold)
  USER_DATA=create_USER_DATA(GS_file)
  
  enrich.go=enricher_internal(gene=sig_gene_symbol,
                              universe = bg_gene_symbol,
                              pvalueCutoff = 1,
                              pAdjustMethod = 'BH',
                              minGSSize = 1,
                              maxGSSize = 1000,
                              qvalueCutoff = 1,
                              USER_DATA = USER_DATA)
  
  ora_df=data.frame(enrich.go)
  ora_df$GeneRatio=paste0("'",ora_df$GeneRatio) #otherwise it will become a date
  write.csv(ora_df,sprintf('%s/DEV_cortex_thres%s_ORA_results_%s.csv',result_path,pSI.threshold,r2test))
}

# 4.4. GSEA/GCEA analysis
source(paste(code_path,'gsea_functions.R',sep = '/'))
result_path=file.path(data_path,'all_results/analysis4-GSEA')
dir.create(result_path)

r2test=0.4

gene_data=read.csv(sprintf('%s/allgenes_stable_r%s.csv',code_path,r2test))%>%rename_with(~paste0('gene_',.),-c('label'))
  
# merge all data by label and remove regions with NA values
df_all=merge(merge(gene_data,pc1,by='label'),lr_centroid,by='label')%>%filter(complete.cases(.))

# do correlation between PC1 and gene expression
gene_df=df_all %>% dplyr::select(starts_with('gene_'))
# get centroid.l/centroid.r for 66 regions
centroid.l=df_all %>% filter(grepl('^L_',label)) %>%  dplyr::select(starts_with('centroid'))
centroid.r=df_all %>% filter(grepl('^R_',label)) %>%  dplyr::select(starts_with('centroid'))
brain_df=df_all %>% dplyr::select(starts_with('brain'))

#perm.id=rotate.parcellation(as.matrix(centroid.l),as.matrix(centroid.r),nrot=10000)
#saveRDS(perm.id,paste(code_path,'66perm_id.rds',sep='/'))
perm.id=readRDS(paste(code_path,'66perm_id.rds',sep='/')) 
  
 null_brain_df=build_null_brain_df(perm.id,brain_df)
  
# get r value for analysis
true_r_vals=cor(gene_df,brain_df,method='pearson')
null_r_vals=cor(gene_df,null_brain_df,method='pearson')
  
  
## DO GSEA/GCEA for GO- CC/MF/BP enrichment analysis
for (ont2use in c('CC','MF','BP')){
  USER_DATA=get_GO_data('org.Hs.eg.db',ont2use,'SYMBOL') # or
  
  gsea_df=enrich_with_snull(true_r_vals,
                            null_r_vals,
                            USER_DATA,
                            minGSSize=10,
                            maxGSSize=200,
                            by='gsea')
  
  gcea_df=enrich_with_snull(true_r_vals,
                            null_r_vals,
                            USER_DATA,
                            minGSSize=10,
                            maxGSSize=200,
                            by='gcea')
  
  write.csv(gsea_df,sprintf('%s/GO_%s_GSEA_results_%s_10-200.csv',result_path,ont2use,r2test))
  write.csv(gcea_df,sprintf('%s/GO_%s_GCEA_results_%s.csv',result_path,ont2use,r2test))
}
  
  
# DO GSEA/GCEA for developmental stages
  for (pSI.threshold in c(0.01)){ # 0.01 is used given the number of genes after intersection with gene expression data
    GS_file=sprintf('%s/pSI/pSI_only_cortex_table_%s.csv',code_path,pSI.threshold)
    USER_DATA=create_USER_DATA(GS_file)
    
    gsea_df=enrich_with_snull(true_r_vals,
                              null_r_vals,
                              USER_DATA,
                              minGSSize=10,
                              maxGSSize=1000,
                              by='gsea')
    
    gcea_df=enrich_with_snull(true_r_vals,
                              null_r_vals,
                              USER_DATA,
                              minGSSize=10,
                              maxGSSize=1000,
                              by='gcea')
    
    write.csv(gsea_df,sprintf('%s/DEV_cortex_%s_GSEA_results_%s.csv',result_path,pSI.threshold,r2test))
    write.csv(gcea_df,sprintf('%s/DEV_cortex_%s_GCEA_results_%s.csv',result_path,pSI.threshold,r2test))
}

# 4.5. Plot significant GO terms
for (ont2use in c('CC','MF','BP')){
  USER_DATA=get_GO_data('org.Hs.eg.db',ont=ont2use,'SYMBOL')
  d = godata('org.Hs.eg.db', ont=ont2use,keytype = "SYMBOL",
             computeIC=FALSE) # computeIC=FALSE when use Wang method
  gcea_res_list=enrich_with_snull(true_r_vals,
                                  null_r_vals,
                                  USER_DATA,
                                  minGSSize=10,
                                  maxGSSize=200,
                                  by='gcea',
                                  returnDF = F) # this will return the obj and res table
  gcea_res_updated=update_gcea_res.obj(gcea_res_list,
                                       settype = ont2use,
                                       filter_fdrp=T,
                                       threshold=0.05)
  gcea_res_updated=pairwise_termsim(gcea_res_updated, method="Wang", semData = d)
  p=treeplot(gcea_res_updated,
             showCategory = 30,
             cex_category = 0.6,
             nCluster = 6,
             nWords = 3,
             color='pvalue',
             offset = rel(14))+
    scale_size(guide = "none")
  ggsave(sprintf('%s/GO_%s_GCEA_plot_%s.png',result_path,ont2use,r2test),
         p, width = 16, height = 8)
}
```

