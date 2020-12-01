# read p values for all the substance
df_list=list()
for (subname in c('alc','nic','coc','met','can')){
  file2read=paste('Model2_results_',subname,'.xls',sep = '')
  df_list[[subname]]=read.csv(file2read,sep = '\t')%>%mutate(drug=rep(subname,77))
}
all_df=bind_rows(df_list)


# adjust p values for CT/SA/subcortial volume across substance types
ct_adjust=all_df %>% filter(grepl("thick|Thickness",region)) %>% mutate(p_adj_cross_substance=p.adjust(p_diag,method="fdr"))
sa_adjust=all_df %>% filter(grepl("surf|SurfArea",region)) %>% mutate(p_adj_cross_substance=p.adjust(p_diag,method="fdr"))
vo_adjust=all_df %>% filter(!grepl("surf|thick|SurfArea|Thickness",region)) %>% mutate(p_adj_cross_substance=p.adjust(p_diag,method="fdr"))
all_df_mcc=rbind(ct_adjust,sa_adjust,vo_adjust)%>%mutate(mcc_cross_substance=case_when(p_adj_cross_substance<0.05~'yes',
                                                                            TRUE~'no'))

# save results
for (subname in c('alc','nic','coc','met','can')){
  savename=paste('Model2_results_',subname,'_new_mcc.csv',sep = '')
  df2save=all_df_mcc%>%filter(drug==subname)
  df2save=df2save[order(match(df2save$region,df_list[[subname]]$region)),] # match the orignal order
  row.names(df2save) <- NULL
  write.csv(df2save,savename)
  rm(df2save)}
