# longitidunal plot
library("plyr")
library("lattice")
library("ggplot2")
library("dplyr")
library("readr")
library("rmarkdown")
library("Rmisc")
library("devtools")
library("gghalves")
library("lmerTest")
library("sjPlot")
library("gridExtra")
library("grid")
longitudinal_plot <- function(d,xlab2show,ylab2show,title2show){
f5 <- ggplot(data = d, aes(y = y)) +
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = 'dodgerblue', size = 1.5,
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = 'darkorange', size = 1.5,alpha = .6) +
  geom_line(aes(x = xj, group = id), color = 'lightgray', alpha = .3) +
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.25),
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
    fill = 'dodgerblue') +
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .15),
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
    fill = 'darkorange') +
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3),
    side = "l", fill = 'dodgerblue') +
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3),
    side = "r", fill = "darkorange") +
  #Define additional settings
  scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0, 3),expand = c(0,0)) +
  xlab(xlab2show) + ylab(ylab2show) +
  ggtitle(title2show) +
  theme_classic()
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max))
f5
return(f5)
}



longitudinal_plot_by_sep <- function(d,xlab2show,ylab2show,title2show){
 f5 <-  ggplot(data=d, aes(y=y)) +
    #Add geom_() objects
    geom_point(data = d %>% filter(sep=='0') %>% filter(x=="1") , aes(x=xj), color = 'dodgerblue', size = 1.5,
               alpha = .6) +
    geom_point(data = d %>% filter(sep=='0') %>% filter(x=="2"), aes(x=xj), color = 'dodgerblue', size = 1.5,
               alpha = .6) +
    geom_point(data = d %>% filter(sep=='1')%>% filter(x=="1") , aes(x=xj), color = 'darkorange', size = 1.5,
               alpha = .6) +
    geom_point(data = d %>% filter(sep=='1')%>% filter(x=="2"), aes(x=xj), color = 'darkorange', size = 1.5,
               alpha = .6) +
    geom_line(aes(x=xj, group=id,colour=sep), alpha = .3) +
   scale_color_manual(values =c('dodgerblue','darkorange'))+
    
    geom_half_boxplot(
      data = d %>% filter(sep=='0')%>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.25),
      side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
      fill = 'dodgerblue',alpha=0.6) +
    geom_half_boxplot(
      data = d %>% filter(sep=='0')%>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .15),
      side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
      fill = 'dodgerblue',alpha=0.6) +
    
    geom_half_violin(
      data = d %>% filter(sep=='0')%>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3),
      side = "l", fill = 'dodgerblue',alpha=0.6) +
    geom_half_violin(
      data = d %>% filter(sep=='0')%>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3),
      side = "r", fill = "dodgerblue",alpha=0.6) +
    geom_half_boxplot(
      data = d %>% filter(sep=='1')%>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.25),
      side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
      fill = 'darkorange',alpha=0.6) +
    geom_half_boxplot(
      data = d %>% filter(sep=='1')%>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .15),
      side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2,
      fill = 'darkorange',alpha=0.6) +
    geom_half_violin(
      data = d %>% filter(sep=='1')%>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3),
      side = "l", fill = 'darkorange',alpha=0.6) +
    geom_half_violin(
      data = d %>% filter(sep=='1')%>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3),
      side = "r", fill = "darkorange",alpha=0.6) +
    
    
    scale_x_continuous(breaks=c(1,2), labels=c("T1", "T2"), limits=c(0, 3),expand = c(0,0)) + # need to add expand to turn off padding
    xlab(xlab2show) + ylab(ylab2show) +
    ggtitle(title2show) +
    theme_classic()
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max))
  f5
return(f5)
  }

longitudinal_plot_by_site <- function(orig_df,col2plot,xlab2show,ylab2show,ylim_min,ylim_max){
plot_list=list()
for (n in c(1:8)){
  df2plot=orig_df %>% filter(all_site==n)
  d=data.frame(x=df2plot$time,
               y=df2plot[,col2plot],
               id=df2plot$subcode,
               sex=df2plot$sex)
  d$xj=jitter(d$x,amount=0.09)
  d$y=jitter(d$y,amount=0.5)
  f=longitudinal_plot(d,xlab2show,ylab2show,sprintf('Site %d',n))+coord_cartesian(ylim=c(ylim_min, ylim_max))
  plot_list[[n]]=f
}
main=grid.arrange(grobs=plot_list,ncol=4)
}



# 
# 
# # Plot MCQ
# df2plot=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/MCQ_kirby/final_MCQ_cov.csv')
# numeric_col=c('age','mod_pds','mean_FD', 'Kest','time')
# factor_col=c('subcode','handedness','sex','all_site')
# 
# 
# for(col in factor_col){
#   df2plot[,col]=as.factor(as.character(df2plot[,col]))
# }
# for (col in numeric_col){
#   df2plot[,col]=as.numeric(as.character(df2plot[,col]))
# }
# col2use=append(numeric_col,factor_col)
# df2plot=df2plot[,col2use]
# 
# 
# test_model=lmer(Kest~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type='pred')$age
# 
# d=data.frame(x=df2plot$time,
#              y=df2plot$Kest,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.05)
# longitudinal_plot(d,'Time','K value','MCQ k value')
# 
# longitudinal_plot_by_sex(d,'Time','K value','MCQ K value by sex (BLue is male)')
# longitudinal_plot_by_site(df2plot,'Kest','Time','Kest')
# 
# 
# 
# # PLOT SSRT
# df2plot=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/SSRT/final_SSRT_cov.csv')
# numeric_col=c('age','mod_pds','mean_FD', 'SSRT','time')
# factor_col=c('subcode','handedness','sex','all_site')
# 
# 
# for(col in factor_col){
#   df2plot[,col]=as.factor(as.character(df2plot[,col]))
# }
# for (col in numeric_col){
#   df2plot[,col]=as.numeric(as.character(df2plot[,col]))
# }
# col2use=append(numeric_col,factor_col)
# df2plot=df2plot[,col2use]
# 
# test_model=lmer(SSRT~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# 
# d=data.frame(x=df2plot$time,
#              y=df2plot$SSRT,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.05)
# longitudinal_plot(d,'Time','SSRT','SSRT')
# 
# longitudinal_plot_by_sex(d,'Time','SSRT','SSRT by sex (BLue is male)')
# longitudinal_plot_by_site(df2plot,'SSRT','Time','SSRT',0,300)
# 
# # PLOT AUDIT BSL012
# df2plot=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT/final_AUDIT_cov_bsl012.csv')
# numeric_col=c('age','mod_pds','mean_FD', 'audit_total','time')
# factor_col=c('subcode','handedness','sex','all_site')
# 
# 
# for(col in factor_col){
#   df2plot[,col]=as.factor(as.character(df2plot[,col]))
# }
# for (col in numeric_col){
#   df2plot[,col]=as.numeric(as.character(df2plot[,col]))
# }
# col2use=append(numeric_col,factor_col)
# df2plot=df2plot[,col2use]
# 
# library("lmerTest")
# test_model=lmer(audit_total~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type = 'pred')$age
# 
# d=data.frame(x=df2plot$time,
#              y=df2plot$audit_total,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.5)
# longitudinal_plot(d,'Time','Audit total','Audit total')
# 
# longitudinal_plot_by_sex(d,'Time','Audit total','Audit total by sex (BLue is male)')
# longitudinal_plot_by_site(df2plot,'audit_total','Time','Audit total',0,30)
# 
# 
# 
# # plot cantab
# df2plot=read.csv('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/intersect_cov/cantab/final_CGT_cov.csv')
# numeric_col=c('age','mod_pds','mean_FD','time',"CGT.Delay.aversion", "CGT.Deliberation.time", "CGT.Overall.proportion.bet","CGT.Quality.of.decision.making")
# factor_col=c('subcode','handedness','sex','all_site')
# 
# 
# for(col in factor_col){
#   df2plot[,col]=as.factor(as.character(df2plot[,col]))
# }
# for (col in numeric_col){
#   df2plot[,col]=as.numeric(as.character(df2plot[,col]))
# }
# col2use=append(numeric_col,factor_col)
# df2plot=df2plot[,col2use]
# 
# 
# test_model=lmer(CGT.Delay.aversion~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type = 'pred')$age
# 
# 
# d=data.frame(x=df2plot$time,
#              y=df2plot$CGT.Delay.aversion,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.5)
# longitudinal_plot(d,'Time','CGT.Delay.aversion','CGT.Delay.aversion')
# 
# longitudinal_plot_by_sex(d,'Time','Audit total','Audit total by sex (BLue is male)')
# 
# 
# 
# x
# 
# # 
# library(readxl)
# df2plot=data.frame(read_excel('/Users/zhipeng//Google Drive//post-doc/MID_network_FU2/process_code/fc_gt_analysis/network_extraction//LME model and plot//all_gt_measures_cov.xlsx'))
# numeric_col=c('age','mod_pds','mean_FD','time',"aCp" ,"aLp", "alocE" , "agE" , "adeg", "abw" ,"amod"   )
# factor_col=c('subcode','handedness','sex','all_site')
#    
# 
# for(col in factor_col){
#   df2plot[,col]=as.factor(as.character(df2plot[,col]))
# }
# for (col in numeric_col){
#   df2plot[,col]=as.numeric(as.character(df2plot[,col]))
# }
# col2use=append(numeric_col,factor_col)
# df2plot=df2plot[,col2use]
# test_model=lmer(agE~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type = 'pred')$age
# 
# d=data.frame(x=df2plot$time,
#              y=df2plot$agE,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.5)
# longitudinal_plot(d,'Time','Global Efficiency','Global Efficiency')
# longitudinal_plot_by_sex(d,'Time','Global Efficiency','Global Efficiency')
# longitudinal_plot_by_site(df2plot,'agE','time','Global Efficiency',-1,2)
# 
# 
# 
# test_model=lmer(amod~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type = 'pred')$age
# d=data.frame(x=df2plot$time,
#              y=df2plot$amod,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.5)
# longitudinal_plot(d,'Time','MOD','MOD')
# longitudinal_plot_by_sex(d,'Time','Global Efficiency','Global Efficiency')
# 
# 
# 
# test_model=lmer(aLp~age+sex+all_site+handedness+(1|subcode),df2plot)
# summary(test_model)
# plot_model(test_model,type = 'pred')$age
# 
# 
# 
# longitudinal_plot_by_site(df2plot,'aLp','time','alp',-0.2,2)
# 
# 
# 
# orig_df=df2plot
# 
# 
# fig_list=list()
# for (n in c(1:8)){
#   print(n)
# df2plot=orig_df %>% filter(all_site==n)
# d=data.frame(x=df2plot$time,
#              y=df2plot$aLp,
#              id=df2plot$subcode,
#              sex=df2plot$sex)
# d$xj=jitter(d$x,amount=0.09)
# d$y=jitter(d$y,amount=0.5)
# f=longitudinal_plot(d,'Time','Lp',str(n))
# fig_list[[n]]=f
# }
# main=grid.arrange(grobs=fig_list,ncol=4)
# longitudinal_plot_by_sex(d,'Time','alp','alp')

