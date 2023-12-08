library(dplyr)
library(stringr)
library(qwraps2)
options(qwraps2_markup = 'markdown')
library(rmarkdown)

# prepare ABCD data
setwd('F:/Google Drive/post-doc/p_factor/ABCD_RSI')
rsi_df=read.csv('rsi_gm_dk_DEAP-data-download.csv')



get_rsi_data <- function(rsi_df, 
                             rsi_name=c('n0','nd','nt','n0s2','nds2','nts2'), 
                             timepoint=c('T1','T2')){
  
  demo_vars=c('src_subject_id','sex_at_birth','event_name','age','abcd_site','race_ethnicity','iqc_dmri_total_passqc',
              'pubertdev_ss_male_category','pubertdev_ss_female_category','smri_vol_subcort.aseg_intracranialvolume',
              'fsqc_qc')
  ordered_variable=read.csv('F:/Google Drive/post-doc/p_factor/clean_code/dk_standard_var_order.csv')[,'label']
  txt2replace=paste('dmri_rsi.',rsi_name,'.gm_cort.desikan_',sep = '')
  # if 2 is not in the measure, need to exclude those end with 2
  if (!(grepl('2',rsi_name,fixed=T))){
  tmp_df=rsi_df %>% select(all_of(demo_vars)|contains(rsi_name)&!contains('2'))
  } else{
    tmp_df=rsi_df %>% select(all_of(demo_vars)|contains(rsi_name))}
  tmp_df=tmp_df %>%rename_if(startsWith(names(.),txt2replace), 
                 ~str_replace(.,txt2replace, ''))%>%
    rename_if(endsWith(names(.),".lh"),  ~str_c('L_',.)) %>%
    rename_if(endsWith(names(.),".rh"),  ~str_c('R_',.))
  colnames(tmp_df) <- gsub('.lh','',colnames(tmp_df),fixed=T) # must add fixed!
  colnames(tmp_df) <- gsub('.rh','',colnames(tmp_df),fixed=T)
  colnames(tmp_df) <- gsub('smri_vol_subcort.aseg_intracranialvolume','ICV',colnames(tmp_df))  
  if (timepoint=='T1'){
      df=tmp_df %>% filter(event_name=='baseline_year_1_arm_1'& iqc_dmri_total_passqc<2&fsqc_qc=='accept')%>%
      mutate(Age=age/12)%>%
      mutate(Sex=case_when(sex_at_birth=='M'~0,
                           sex_at_birth=='F'~1))%>%
      mutate(male_pds=pubertdev_ss_male_category-2)%>% #adjust by the mod (2 is the mod for male)
      mutate(female_pds=pubertdev_ss_female_category-3)%>% #adjust by the mod (3 is the mode for female)
      mutate(mod_PDS=coalesce(male_pds,female_pds))%>% 
      mutate(BSL_pds=coalesce(pubertdev_ss_male_category,pubertdev_ss_female_category)) %>% 
      rename(Site=abcd_site)%>%
      rename(Ethnicity=race_ethnicity)%>%
      rename(SubjID=src_subject_id)%>%
      select(c(SubjID, Age, Sex, Site, Ethnicity,mod_PDS, BSL_pds, ICV)|starts_with('L_')|starts_with('R_'))%>%
      select(-ends_with('_mean'))%>%
      filter(complete.cases(.))
      #reoder the columns so that they are matched with others [the standard FS output order]
      df_demo=df%>%select(c(SubjID, Age, Sex, Site, Ethnicity,mod_PDS, BSL_pds, ICV))
      df_data=df[,as.character(ordered_variable)]
      df.final=cbind(df_demo,df_data)
  } else if (timepoint=='T2'){
    df=tmp_df %>% filter(event_name=='2_year_follow_up_y_arm_1'& iqc_dmri_total_passqc<2&fsqc_qc=='accept')%>%
      mutate(Age=age/12)%>%
      mutate(Sex=case_when(sex_at_birth=='M'~0,
                           sex_at_birth=='F'~1))%>%
      rename(Site=abcd_site)%>%
      rename(Ethnicity=race_ethnicity)%>%
      rename(SubjID=src_subject_id)%>%
      select(c(SubjID, Age, Sex, Site, Ethnicity, ICV)|starts_with('L_')|starts_with('R_'))%>%
      select(-ends_with('_mean'))%>%
      filter(complete.cases(.))
    
      df_demo=df%>%select(c(SubjID, Age, Sex, Site, Ethnicity, ICV))
      df_data=df[,as.character(ordered_variable)]
      df.final=cbind(df_demo,df_data)
  }
  return(df.final)
}

# split data by measures and time point
rsi_measures=c('n0','nd','nt','n0s2','nds2','nts2')
timepoints=c('T1','T2')

for (measure_i in c(1:length(rsi_measures))){
  for (time_i in c(1:length(timepoints))){
    rsi_name=rsi_measures[measure_i]
    time_name=timepoints[time_i]
    save_name=sprintf('ABCD_%s_%s_ready.csv',time_name, rsi_name)
    data=get_rsi_data(rsi_df,rsi_name,time_name)
    write.csv(data,save_name, row.names = F)
  }
}


