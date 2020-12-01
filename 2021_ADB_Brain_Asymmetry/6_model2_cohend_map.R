library(ggplot2)
library(ggseg)
library(gridExtra)

# load list of regions
area_names=read.csv('regions.csv',header = FALSE)

# create table to plot
table2plot <- data.frame(matrix(nrow=length(colnames(aidata)[grep("AI", colnames(aidata))]), ncol=0),stringsAsFactors = F)
table2plot[,"region"] <- sub("AI_", "", colnames(aidata)[grep("AI", colnames(aidata))])
table2plot[,"cohensd"] <- round(as.numeric(d_dx),digits=3)
table2plot[,"area"]=area_names$V3


# make aseg plot
aseg_data=table2plot[-grep("surf|thick|Thickness|SurfArea", table2plot$region),]
aseg_data$area=as.character(aseg_data$area)
p1=ggseg(.data=aseg_data,atlas=aseg,mapping=aes(fill=cohensd), colour="darkgrey",hemisphere='left',size=.8,alpha=0.8)+
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ,limits=c(-0.2,0.2),aes(title="Conhen's d"))+
  ggtitle('C) Subcortical Volume')+
  theme(axis.title = element_blank(),
        axis.text = element_blank())


#make surfarea plot
surface_data=table2plot[grep("surfavg",table2plot$region),] #total SA was excluded by the list in regions.csv
surface_data$area=as.character(surface_data$area)
p2=ggseg(.data=surface_data,mapping=aes(fill=cohensd),colour="darkgrey",hemisphere='left',alpha=0.8)+
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ,limits=c(-0.2,0.2),aes(title="d"))+
  ggtitle('B) Surface Area')+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')


#make thickness plot
thickness_data=table2plot[grep("thickavg",table2plot$region),]
thickness_data$area=as.character(thickness_data$area)
p3=ggseg(.data=thickness_data,mapping=aes(fill=cohensd),colour="darkgrey",hemisphere='left',alpha=0.8)+
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ,limits=c(-0.2,0.2),aes(title="d"))+
  ggtitle('A) Cortical Thickness')+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none')

p4 <- grid.arrange(p3,p2,p1, widths=c(2,2,0.5),layout_matrix=rbind(c(1,3),
                                                                   c(1,3),
                                                                   c(2,3),
                                                                   c(2,3)))
str2save=substring(substance2test,nchar(substance2test)-2)
ggsave(paste('Model2_cohend_map_',str2save,'.pdf',sep = ''),p4)


