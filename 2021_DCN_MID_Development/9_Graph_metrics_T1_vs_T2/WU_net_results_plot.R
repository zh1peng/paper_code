library(gridExtra)
library(grid)
library(dplyr)
library(lmerTest)
library(sjPlot)
library(ggplot2)
setwd('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/net plot')
source('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/longitudinal_plot.R')
df=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_analysis_new_subid/net plot//data2line_plot_mean_se.csv')
df$time=as.factor(df$time)
df$sparsity=5+2*(df$sparsity-1)

gt_measures=read.csv('WU_all_gt_measures.csv')
bsl_good_behav=read.csv('bsl_RT_good_1491indx.csv')
fu2_good_behav=read.csv('fu2_RT_good_1365indx.csv')
good_idx_df=rbind(bsl_good_behav,fu2_good_behav)
gt_measures=cbind(gt_measures,good_idx_df) %>% filter(good_behav==1)
node_info=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_extraction_new_bu//Nodal analysis//nodes166_info.csv')
bsl_cov=subset(read.csv('bsl_cov1304.csv'),select=-c(X,good_behav))
fu2_cov=subset(read.csv('fu2_cov1240.csv'),select=-c(X,good_behav))
cov_df=rbind(bsl_cov,fu2_cov)
# colnames(cov_df)



gt_df=cbind(cov_df,gt_measures)
numeric_col=c('age','mod_pds','mean_FD', 'time','age_diff') # time is put as numeric_col for plot purpose. but in the analysis it was included as factor!
numeric_col=c(numeric_col,colnames(gt_df)[grep("node_|net_",colnames(gt_df))])
factor_col=c('subcode','handedness','sex','all_site')

for(col in factor_col){
  gt_df[,col]=as.factor(as.character(gt_df[,col]))
}
for (col in numeric_col){
  gt_df[,col]=as.numeric(as.character(gt_df[,col]))
}
gt_df$time_f=as.factor(gt_df$time) 


# plot Lp
p11=ggplot(df, aes(x=sparsity, y=net_Lp_mean, group=time, color=time)) + 
  geom_errorbar(aes(ymin=net_Lp_mean-net_Lp_se, ymax=net_Lp_mean+net_Lp_se), width=.2) +
  geom_line() + geom_point()+
  labs(x="Sparsity %", y = "Shortest Path Length")+
  xlim(1,40)+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("T1", "T2"),
                     name='Time')

d=data.frame(x=gt_df$time,
             y=gt_df$net_aLp,
             id=gt_df$subcode,
             sep=gt_df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p12=longitudinal_plot(d,'Time','AUC','')+theme_classic()


lme_model=lmer(net_aLp~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
p13=plot_model(lme_model,type = 'pred',axis.title = c('Time','Marginal Effects'), axis.labels = c(c('BSL','FU2')),
           title = '')$time_f+scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0.5, 2.5))+
            theme_classic()


longitudinal_plot_by_sep(d,'Time','AUC','')+theme_classic()
lme_model=lmer(net_aLp~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
plot_model(lme_model,type = 'pred')$sex
summary(lme_model)

# plot Degree
p21=ggplot(df, aes(x=sparsity, y=net_deg_mean, group=time, color=time)) + 
  geom_errorbar(aes(ymin=net_deg_mean-net_deg_se, ymax=net_deg_mean+net_deg_se), width=.2) +
  geom_line() + geom_point()+
  labs(x="Sparsity %", y = "Strength")+
  xlim(1,40)+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("T1", "T2"),
                     name='Time')

d=data.frame(x=gt_df$time,
             y=gt_df$net_adeg,
             id=gt_df$subcode,
             sex=gt_df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p22=longitudinal_plot(d,'Time','AUC','')+theme_classic()


lme_model=lmer(net_adeg~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
p23=plot_model(lme_model,type = 'pred',title='', axis.title = c('Time','Marginal Effects'))$time_f+
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0.5, 2.5))+
  theme_classic()

# plot Lp
p31=ggplot(df, aes(x=sparsity, y=net_gE_mean, group=time, color=time)) + 
  geom_errorbar(aes(ymin=net_gE_mean-net_gE_se, ymax=net_gE_mean+net_gE_se), width=.2) +
  geom_line() + geom_point()+
  labs(x="Sparsity %", y = "Global Efficiency")+
  xlim(1,40)+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("T1", "T2"),
                     name='Time')

d=data.frame(x=gt_df$time,
             y=gt_df$net_agE,
             id=gt_df$subcode,
             sex=gt_df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p32=longitudinal_plot(d,'Time','AUC','')+theme_classic()
lme_model=lmer(net_agE~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
p33=plot_model(lme_model,type = 'pred',title='', axis.title = c('Time','Marginal Effects'))$time_f+
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0.5, 2.5))+
  theme_classic()



# plot Cp
p41=ggplot(df, aes(x=sparsity, y=net_Cp_mean, group=time, color=time)) + 
  geom_errorbar(aes(ymin=net_Cp_mean-net_Cp_se, ymax=net_Cp_mean+net_Cp_se), width=.2) +
  geom_line() + geom_point()+
  labs(x="Sparsity %", y = "Clustering Coefficient")+
  xlim(1,40)+
  theme_classic()+
  scale_color_manual(values=c('dodgerblue','darkorange'),
                     labels = c("T1", "T2"),
                     name='Time')

d=data.frame(x=gt_df$time,
             y=gt_df$net_aCp,
             id=gt_df$subcode,
             sex=gt_df$sex)
d$xj=jitter(d$x,amount=0.09)
d$y=jitter(d$y,amount=0.05)
p42=longitudinal_plot(d,'Time','AUC','')+theme_classic()
lme_model=lmer(net_aCp~time_f+sex+all_site+handedness+mean_FD+mod_pds+age_diff+(1|subcode),data=gt_df)
p43=plot_model(lme_model,type = 'pred',title='', axis.title = c('Time','Marginal Effects'))$time_f+
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0.5, 2.5))+
  theme_classic()


main= grid.arrange(p21,p22,p23,p11,p12,p13,p41,p42,p43,ncol=3)
ggsave("net_results.tiff", main,width = 12,height = 7)





