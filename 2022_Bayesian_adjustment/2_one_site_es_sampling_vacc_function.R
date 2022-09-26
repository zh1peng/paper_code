#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
idx=as.numeric(args[1])
source('/gpfs1/home/z/c/zcao4/Bayes_analysis/analysis_functions.R')
data_df=read.csv('/gpfs1/home/z/c/zcao4/Bayes_analysis/combat_data.csv')


numeric_col=colnames(data_df)[grep('L_|R_|Age|ICV',colnames(data_df))]
factor_col=c('Sex','Dependentanydrug','Site')


for(col in factor_col){
  data_df[,col]=as.factor(as.character(data_df[,col]))
}
for (col in numeric_col){
  data_df[,col]=as.numeric(as.character(data_df[,col]))
}


save_dir='/gpfs1/home/z/c/zcao4/Bayes_analysis/site_data'
site=unique(data_df$Site)
site_i=site[idx] 
save_name1=sprintf('%s/site_%s_ob_es.csv',save_dir,site_i)
save_name2=sprintf('%s/site_%s_sampled_es.csv',save_dir,site_i)
site_df=data_df %>% filter(Site==site_i) %>% droplevels

site_es=get_site_es(site_df)
site_sampled_es=get_sampled_es_for_site(data_df,site_df,1000)
write.csv(site_es,save_name1)
write.csv(site_sampled_es,save_name2)





