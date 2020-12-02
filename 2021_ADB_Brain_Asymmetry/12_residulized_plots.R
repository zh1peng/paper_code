library(ggplot2)
library(ggpubr)
library(plyr)
library(gridExtra)


# left and right plot
df2plot=data.frame(NAcc=c(as.numeric(resid_LRdata$L_accumb),as.numeric(resid_LRdata$R_accumb)),
                   hem=c(rep('Left',dim(resid_LRdata)[1]),rep('Right',dim(resid_LRdata)[1])),
                   anydrug=rbind(resid_LRdata,resid_LRdata)$Dependentanydrug,
                   drugtype=rbind(resid_LRdata,resid_LRdata)$PrimaryDrug)



df2plot$ifdependent=mapvalues(df2plot$anydrug, from=c('0','1'),to=c('Nondependent','Dependent'))
 

##===================  Model 1 plot ============================================ 
  p1=ggplot(data = df2plot,aes(x = ifdependent, y = NAcc, fill = hem))+
  scale_fill_manual("legend", values = c("Left" = "#1d91c0", "Right" = "#c7e9b4"))+
  geom_point(shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
    ylab('Residulized NAcc Volume')

resid_aidata$ifdependent=mapvalues(resid_aidata$Dependentanydrug, from=c('0','1'),to=c('Nondependent','Dependent'))
  p2=ggplot(data = resid_aidata,aes(x = ifdependent, y = AI_accumb,fill="AI"))+
    scale_fill_manual("legend", values = c("AI" = "#253494"))+
    geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
    geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  theme(legend.position = "none",
        axis.title.x=element_blank())+
    ylab('Residulized NAcc AI')
   
  
  pdf('Model1_residulized_plot.pdf') 
  grid.arrange(p1,p2,ncol=1) 
  dev.off()


    
    
##===================  Model 2 plot ============================================
df2plot$substance=mapvalues(df2plot$drugtype,from=c('0','1','2','3','4','5'), 
                            to=c('Nondependent','Alcohol','Nicotine','Cocaine','Methamphetamine','Cannabis'))
p3=  ggplot(data = df2plot,aes(x = substance, y = NAcc, fill = hem))+
    scale_fill_manual("legend", values = c("Left" = "#1d91c0", "Right" = "#c7e9b4"))+
    geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
    geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank())+
      ylab('Residulized NAcc Volume')
  
  
  resid_aidata$substance=mapvalues(resid_aidata$PrimaryDrug,from=c('0','1','2','3','4','5'), 
                                   to=c('Nondependent','Alcohol','Nicotine','Cocaine','Methamphetamine','Cannabis'))

 p4= ggplot(data = resid_aidata,aes(x = substance,  y = AI_accumb,fill="AI"))+
    scale_fill_manual("legend", values = c("AI" = "#253494"))+
    geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
    geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
    geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x = element_text(size=8.5,colour = 'black'))+
    ylab('Residulized NAcc AI')
 
 pdf('Model2_residulized_plot.pdf') 
 grid.arrange(p3,p4,ncol=1) 
  dev.off()

