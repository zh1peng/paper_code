library(parameters)
library(dplyr)

cohen_es <- function(df, t, n1, n2){
  d <- t*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
  return(d)}

get_site_es <- function(site_df){
site_es=list()
for (col in colnames(site_df)[grep('L_|R_',colnames(site_df))]){
  f2use =sprintf('%s~Dependentanydrug+Sex+Age+ICV',col)
  if (nlevels(site_df$Sex)==1){
    f2use =sprintf('%s~Dependentanydrug+Age+ICV',col)
  }
  site_lm=lm(as.formula(f2use), data=site_df)
  g1_n=sum(site_df$Dependentanydrug==0)
  g2_n=sum(site_df$Dependentanydrug==1)
  t_value=parameters(site_lm)[2,6]
  df_value=parameters(site_lm)[2,7]
  site_es[[col]]=cohen_es(df_value,t_value,g1_n,g2_n)
}
return(data.frame(site_es))
}


#=============================calculate sampled es==============================
get_sampled_es_for_site <- function(data_df,site_df,n_iter){
  all_sampled_es=data.frame()
  for (sample_i in c(1:n_iter)){
  case_m_n=dim(site_df[site_df$Dependentanydrug==1 & site_df$Sex ==1,])[1]
  case_f_n=dim(site_df[site_df$Dependentanydrug==1 & site_df$Sex ==2,])[1]
  hc_m_n=dim(site_df[site_df$Dependentanydrug==0 & site_df$Sex ==1,])[1]
  hc_f_n=dim(site_df[site_df$Dependentanydrug==0 & site_df$Sex ==2,])[1]
  case_m_df=data.frame()
  case_f_df=data.frame()
  hc_m_df=data.frame()
  hc_f_df=data.frame()
  if (case_m_n>0){
    case_m_df=data_df %>% filter (Dependentanydrug==1 & Sex ==1) %>% sample_n(case_m_n)}
  if (case_f_n>0){
    case_f_df=data_df %>% filter (Dependentanydrug==1 & Sex ==2) %>% sample_n(case_f_n)}
  if (hc_m_n>0){
    hc_m_df=data_df %>% filter (Dependentanydrug==0 & Sex ==1) %>% sample_n(hc_m_n)}
  if (hc_f_n>0){
    hc_f_df=data_df %>% filter (Dependentanydrug==0 & Sex ==2) %>% sample_n(hc_f_n)}
  
  sampled_df=rbind(case_m_df,case_f_df,hc_m_df,hc_f_df) %>% droplevels
  sampled_es=get_site_es(sampled_df)
  all_sampled_es=rbind(all_sampled_es,sampled_es)
  }
  return(all_sampled_es)
}



get_sampled_es_for_sample_n <- function(data_df,case_m_n,case_f_n,hc_m_n,hc_f_n,n_iter){
  all_sampled_es=data.frame()
  for (sample_i in c(1:n_iter)){
    case_m_df=data.frame()
    case_f_df=data.frame()
    hc_m_df=data.frame()
    hc_f_df=data.frame()
    if (case_m_n>0){
      case_m_df=data_df %>% filter (Dependentanydrug==1 & Sex ==1) %>% sample_n(case_m_n)}
    if (case_f_n>0){
      case_f_df=data_df %>% filter (Dependentanydrug==1 & Sex ==2) %>% sample_n(case_f_n)}
    if (hc_m_n>0){
      hc_m_df=data_df %>% filter (Dependentanydrug==0 & Sex ==1) %>% sample_n(hc_m_n)}
    if (hc_f_n>0){
      hc_f_df=data_df %>% filter (Dependentanydrug==0 & Sex ==2) %>% sample_n(hc_f_n)}
    sampled_df=rbind(case_m_df,case_f_df,hc_m_df,hc_f_df) %>% droplevels
    sampled_es=get_site_es(sampled_df)
    all_sampled_es=rbind(all_sampled_es,sampled_es)
  }
  return(all_sampled_es)
}

