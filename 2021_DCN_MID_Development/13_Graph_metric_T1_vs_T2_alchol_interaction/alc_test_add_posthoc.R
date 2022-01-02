library(lmerTest)
library(dplyr)
# test with fu2 alc group
df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov//AUDIT analysis new subid//audit_brain_cov.csv')
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/nodes166_info_updated.csv')
numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(df)[grep("bu_|wu_|behav_|audit_|mtg_",colnames(df))])
factor_col=c('subcode','handedness','sex','all_site')

for(col in factor_col){
  df[,col]=as.factor(as.character(df[,col]))
}
for (col in numeric_col){
  df[,col]=as.numeric(as.character(df[,col]))
}
df$time_f=as.factor(df$time) 


# add audit group
bsl_df=df %>% filter(time==1)
fu2_df=df %>% filter(time==2)
fu2_df=fu2_df %>% mutate(fu2_alc=case_when(audit_total>7~1,
                                           audit_total<8~0))
fu2_df$fu2_alc=as.factor(fu2_df$fu2_alc)
bsl_df$fu2_alc=as.factor(fu2_df$fu2_alc)
df=rbind(bsl_df,fu2_df)


time_p=list()
time_t=list()
int_t=list()
int_p=list()
var_t=list()
var_p=list()
int_null_t=list()
var_null_t=list()
for (col1 in c('mgt_vsr_node_MTG')){
  for (col2 in c('fu2_alc')){
    col2show=sprintf('%s-%s',col1,col2)
    tmp_data=df
    f2use=sprintf('%s~mod_pds+handedness+time_f*%s+sex+all_site+age_diff+mean_FD+(1|subcode)',col1,col2)
    lme_model=lmer(as.formula(f2use), data=tmp_data)
    time_t[[col2show]]=summary(lme_model)$coefficients['time_f2','t value']
    time_p[[col2show]]=summary(lme_model)$coefficients['time_f2','Pr(>|t|)']
    var_t[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'t value']
    var_p[[col2show]]=summary(lme_model)$coefficients[sprintf('%s1',col2),'Pr(>|t|)']
    int_t[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'t value']
    int_p[[col2show]]=summary(lme_model)$coefficients[sprintf('time_f2:%s1',col2),'Pr(>|t|)']
    perm_results=node_perm_parallel_function(df,col1,14)
    int_null_t[[col2show]]=perm_results$null_int_t
    var_null_t[[col2show]]=perm_results$null_alc_t
  }
}
net_df=data.frame(time_t=unlist(time_t),
                  time_p=unlist(time_p),
                  var_t=unlist(var_t),
                  var_p=unlist(var_p),
                  int_t=unlist(int_t),
                  int_p=unlist(int_p))
net_df=net_df %>% mutate(fdr_int_p=p.adjust(int_p, method="fdr"))
net_df=net_df %>% mutate(fdr_alc_p=p.adjust(var_p, method="fdr"))
net_df=net_df %>% mutate(fdr_time_p=p.adjust(time_p, method="fdr"))


tmp_quantile=data.frame(l025=unlist(sapply(var_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(var_null_t, function(x){quantile(x,0.975)})))
net_df$alc_l025=tmp_quantile$l025
net_df$alc_u975=tmp_quantile$u975

tmp_quantile=data.frame(l025=unlist(sapply(int_null_t, function(x){quantile(x,0.025)})),
                        u975=unlist(sapply(int_null_t, function(x){quantile(x,0.975)})))
net_df$int_l025=tmp_quantile$l025
net_df$int_u975=tmp_quantile$u975

net_df=net_df %>% mutate(alc_ifsig=(var_t<alc_l025 | var_t>alc_u975))
net_df=net_df %>% mutate(alc_sig_dir=case_when(alc_ifsig==1 & var_t>0 ~1,
                                               alc_ifsig==1 & var_t<0 ~2,
                                               alc_ifsig==0 ~0))
net_df=net_df %>% mutate(int_ifsig=(int_t<int_l025 | int_t>int_u975))
net_df=net_df %>% mutate(int_sig_dir=case_when(int_ifsig==1 & int_t>0 ~1,
                                               int_ifsig==1 & int_t<0 ~2,
                                               int_ifsig==0 ~0))

net_df=cbind(net_df,net_info)
rownames(net_df)=c('mgt_vsr_node_MTG')
write.csv(net_df,'/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/MTG_Interaction plot/alc_net_comp_results.csv')

# make interaction plot
library(interactions)
lme_model=lmer(mgt_vsr_node_MTG~time_f*fu2_alc+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
p1=cat_plot(lme_model, pred = time_f, modx=fu2_alc,geom = "line",
               modx.labels = c("LA", "HA"),legend.main = 'Group',
               pred.labels =c("T1", "T2"),
               x.label="Time",
               y.label="Marginal Effects",
               main.title="Left middle temporal gyrus")+
  theme_classic()

plot_model(lme_model,type='int',title='Time x Alcohol Group',axis.title = c('Time','Marginal Effects'))+
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0, 3))+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("LA", "HA"),
                     name='Group')



d=data.frame(x=df$time,
             y=df$mgt_vsr_node_MTG,
             id=df$subcode,
             sep=df$fu2_alc)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
longitudinal_plot_by_sep(d,'Time','AUC','')+theme_classic()



lme_model=lmer(wu_net_aLp~time_f*fu2_alc+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
p2=cat_plot(lme_model, pred = time_f, modx=fu2_alc,geom = "line",
            modx.labels = c("LA", "HA"),legend.main = 'Group',
            pred.labels =c("T1", "T2"),
            x.label="Time",
            y.label="Marginal Effects",
            main.title="Shortest path length")+
  theme_classic()
  



p32=plot_model(lme_model,type = 'int',title='Time x Alcohol Group',axis.title = c('Time','Marginal Effects'))+
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0, 3))+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("LA", "HA"),
                     name='Group')

main=grid.arrange(p1, p2, nrow = 1)
ggsave('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/MTG_Interaction plot/interaction_plot.tiff',
       main,width = 8,height = 4)



# add post-hocs analysis
library(lsmeans)
library(emmeans)

setwd('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/variability_add_MTG/MTG_Interaction plot')
lme_model=lmer(mgt_vsr_node_MTG~time_f*fu2_alc+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
mgt_vsr_post_hoc=emmeans(lme_model,list(pairwise~time_f*fu2_alc),adjust='fdr')
write.csv(as.data.frame(mgt_vsr_post_hoc$`emmeans of time_f, fu2_alc`),'mgt_post_hoc1.csv')
write.csv(as.data.frame(mgt_vsr_post_hoc$`pairwise differences of time_f, fu2_alc`),'mgt_post_hoc2.csv')

lme_model=lmer(wu_net_aLp~time_f*fu2_alc+mod_pds+handedness+mean_FD+sex+all_site+age_diff+(1|subcode),df)
net_path_post_hoc=emmeans(lme_model,list(pairwise~time_f*fu2_alc),adjust='fdr')
write.csv(as.data.frame(net_path_post_hoc$`emmeans of time_f, fu2_alc`),'net_path_post_hoc1.csv')
write.csv(as.data.frame(net_path_post_hoc$`pairwise differences of time_f, fu2_alc`),'net_path_post_hoc2.csv')





##====================== build null models ========================


shuffle_alc <- function(cov){
  shuffled_cov=cov
  bsl_alc=cov%>%filter(time==1)%>%select(fu2_alc)
  shuffled_alc=bsl_alc[sample(1:nrow(bsl_alc)),]
  shuffled_cov$shuffled_alc=unlist(list(shuffled_alc,shuffled_alc)) # way to concat two factor
  return(shuffled_cov)
}


node_perm_single_fucntion <- function(gt_df,col){
  gt_df_shuffle_alc=shuffle_alc(gt_df)
  y=gt_df[,col]
  null_model=lmer(y~time_f*shuffled_alc+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df_shuffle_alc)
  null_t=summary(null_model)$coefficients['shuffled_alc1','t value']
  null_int=summary(null_model)$coefficients['time_f2:shuffled_alc1','t value']
  return(cbind(null_t,null_int))
} 
library(doParallel)
node_perm_parallel_function <- function(gt_df,col,ncore2use){
  registerDoParallel(cores=ncore2use)
  tmp_results=foreach (perm=1:1000) %dopar% node_perm_single_fucntion(gt_df,col)
  null_alc_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_t'])))
  null_int_t=unlist(lapply(tmp_results, function(i) as.numeric(i[,'null_int'])))
  result_list=list('null_alc_t'=null_alc_t,
                   'null_int_t'=null_int_t)
  return(result_list)
}
