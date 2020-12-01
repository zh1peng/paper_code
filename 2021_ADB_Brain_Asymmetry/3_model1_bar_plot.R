library(ggeffects)
library(ggplot2)
library(gridExtra)
library(nlme)
library(grid)
N.cases.LR=list()
N.controls.LR=list()
########################## Model 1 ###################################
cov = c("Age", "Sex", "Dependentanydrug", "Site")
regioni = names(sig)
region2plot=sub('AI_','',regioni)
region=paste('L_',region2plot,sep = '')
LRdata.complete <- LRdata[,c(cov,region)][which(complete.cases(LRdata[,c(cov,region)])),]
N.cases.LR[[region]] <- length(unique(subset(LRdata.complete, Dependentanydrug =="1")$Subject))
N.controls.LR[[region]] <- length(unique(subset(LRdata.complete, Dependentanydrug =="0")$Subject))
outcome <- LRdata.complete[,region]
L_lme=lme(outcome ~ Dependentanydrug + Sex + Age, random=~1|Site, data = LRdata.complete, method="ML", na.action="na.fail")
L_pred=ggpredict(L_lme, terms = "Dependentanydrug")


region=paste('R_',region2plot,sep = '')
LRdata.complete <- LRdata[,c(cov,region)][which(complete.cases(LRdata[,c(cov,region)])),]
outcome <- LRdata.complete[,region]
R_lme=lme(outcome ~ Dependentanydrug + Sex + Age, random=~1|Site, data = LRdata.complete, method="ML", na.action="na.fail")
R_pred=ggpredict(R_lme, terms = "Dependentanydrug")


region=paste('AI_',region2plot,sep = '')
aidata.complete <- aidata[,c(cov,region)][which(complete.cases(aidata[,c(cov,region)])),]
outcome <- aidata.complete[,region]
ai_lme=lme(outcome ~ Dependentanydrug + Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
ai_pred=ggpredict(ai_lme, terms = "Dependentanydrug")




#library(sjPlot)
#plot_model(ai_lme,type='pred')$Dependentanydrug


df_ctl=data.frame(
  group=c('Nondependent'),
  hem=c('NAcc_L','NAcc_R','NAcc_AI'),
  mean_val=c(L_pred[L_pred$x==0,]$predicted, R_pred[R_pred$x==0,]$predicted,ai_pred[ai_pred$x==0,]$predicted),
  std_val=c(L_pred[L_pred$x==0,]$std.error,R_pred[R_pred$x==0,]$std.error,ai_pred[ai_pred$x==0,]$std.error))


df_su=data.frame(
  group=c('Dependent'),
  hem=c('NAcc_L','NAcc_R','NAcc_AI'),
  mean_val=c(L_pred[L_pred$x==1,]$predicted, R_pred[R_pred$x==1,]$predicted,ai_pred[ai_pred$x==1,]$predicted),
  std_val=c(L_pred[L_pred$x==1,]$std.error,R_pred[R_pred$x==1,]$std.error,ai_pred[ai_pred$x==1,]$std.error))

final_df=rbind(df_ctl,df_su)
df2plot1=final_df[final_df$hem!='NAcc_AI',]

p1=ggplot(df2plot1, aes(fill=hem, y=mean_val, x=group)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual("legend", values = c("NAcc_R" = "#1d91c0", "NAcc_L" = "#c7e9b4"))+
  geom_errorbar(aes(ymin=mean_val-std_val, ymax=mean_val+std_val),
                width=.2,                   
                position=position_dodge(.9))+
  ylab('Estimated Margins (Volumne)')+
  scale_y_continuous( limits = c(0, 700),breaks =c(0,200,400,600))+
  
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

df2plot2=final_df[final_df$hem=='NAcc_AI',]

p2=ggplot(df2plot2, aes( y=mean_val, x=group)) + 
  geom_bar(position="dodge", stat="identity",fill="#253494")+
  geom_errorbar(aes(ymin=mean_val-std_val, ymax=mean_val+std_val),
                width=.2,                   
                position=position_dodge(.9))+
  
  ylab('Estimated Margins (AI)')+
  scale_y_continuous( limits = c(-0.05, 0.01),breaks =c(0.01,0,-0.01,-0.02,-0.03,-0.04,-0.05))+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_text(size=8, face="bold",colour = 'black'))

pdf('Model1_bar_plot.pdf')
grid.arrange(p1,p2, nrow=2)
dev.off()
