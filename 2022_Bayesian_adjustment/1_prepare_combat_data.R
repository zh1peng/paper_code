working_dir='F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model\\'
out <- 'yes'
remove_dataset='no'
txt2extract=""
source(paste(working_dir,"1_load_data.R", sep=""))            ## loads data needed for main analysis ** redundant columns will be removed **


#==============================remove control/case only site==============
# control/case only site
tmp = LR2use %>% group_by(Dependentanydrug, Site) %>% count(Dependentanydrug)
if( length(names(which(table(tmp$Site)==1))) != 0) ## if 1 instead of 2, only cases or only control data is present
{
  case_control_only = names(which(table(tmp$Site)==1))
  print(case_control_only)
}

site2exclude=case_control_only
train_df=LR2use %>% filter(!Site %in% site2exclude) 

test=LR2use%>% filter(Site %in% site2exclude) %>% 
                group_by(Site) %>% 
                distinct(PI)
#==============================report sample characteristics==============

table1=train_df %>%  group_by(Site, Dependentanydrug ) %>% summarise(n=n(),
                                                              mean_age=round(mean(Age),2),
                                                              sd_age=round(sd(Age),2),
                                                              Male=round(sum(Sex==1)/n()*100,2))

drug_dict=train_df %>% group_by(Site) %>%  
                      distinct(PrimaryDrug) %>% 
                      filter(!PrimaryDrug==0) 



PI_dict=train_df %>% group_by(Site) %>%  
                    distinct(PI) %>% 
                    filter(!PI=='')



site_n_df=train_df %>% group_by(Site) %>% summarise(n=n()) 
replace_dict=site_n_df %>% arrange(n)%>% mutate(site_idx=c(1:21))
table1$site_idx=plyr::mapvalues(as.character(table1$Site),
                from = as.character(replace_dict$Site),
                to = replace_dict$site_idx)
table1$PrimaryDrug=plyr::mapvalues(as.character(table1$Site),
                            from = as.character(drug_dict$Site),
                            to = as.character(drug_dict$PrimaryDrug)) # 7,23 have two primary drug

table1$PI=plyr::mapvalues(as.character(table1$Site),
                            from = as.character(PI_dict$Site),
                            to = as.character(PI_dict$PI))

table1=table1 %>% mutate(`n (male%)`=paste0(n,'(',Male,'%)')) %>% 
                  mutate(`Age mean(sd)`=paste0(mean_age,'(',sd_age,')')) %>% 
                  mutate(PrimaryDrugName=recode(PrimaryDrug,
                                                    '1'='ALC',
                                                    '2'='NIC',
                                                    '3'='COC',
                                                    '4'='MET',
                                                    '5'='CAN')) %>% 
                  mutate(Dependentanydrug=recode(Dependentanydrug,
                                                 '0'='Control',
                                                 '1'='Case'))
write.csv(table1,'F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model\\table1.csv',row.names = F)

train_df %>% group_by(Dependentanydrug) %>% summarise(n=n())

#==============================apply combat================================
source(paste(working_dir,"combat_for_ENIGMA_sMRI.R", sep="")) 
train_pheno=train_df[,c('Age','Sex','Site','Dependentanydrug')]
train_dat=train_df[,grep('L_|R_|ICV',colnames(train_df))]
train_modcombat=model.matrix(~as.factor(Dependentanydrug)+as.factor(Sex)+as.numeric(Age),data=train_pheno)
train_batch=droplevels(train_df[,'Site'])  # drop level 
combat_model=combat_fit(train_dat, train_batch, mod = train_modcombat, eb = TRUE, verbose = TRUE)
train_combat_results=combat_apply(combat_model, train_dat, train_batch, mod = train_modcombat, verbose = TRUE)
train_combat_df=cbind(train_combat_results$dat.combat, train_pheno)

write.csv(train_combat_df,'F:\\Google Drive\\post-doc\\Bayesian_Project\\vacc_analysis\\combat_data.csv')

#=============================calculate site es for all ROI===================================

