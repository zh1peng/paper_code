library(ggseg)
library(viridis)
library(ggplot2)
# make illustration plots

# make density function
# r2z <- function(r){
# z = .5 * log((1+r)/(1-r))
# }
# cor2plot=as.data.frame(cor2plot) %>% mutate(r.abs=abs(x),
#                              z=r2z(x))

# plot with base 
# d_res=density(cor2plot$x)
# plot(d_res)
# abline(v=mean(res_r),col="black",lwd=3,lty=2) ## fixed

plot_cor <- function(cor2plot,color2use=c('#FDE725FF','#33638DFF')){
 if(!is.data.frame(cor2plot)){cor2plot=as.data.frame(cor2plot)} 
p=ggplot(cor2plot, aes(x=x)) +
  stat_density(geom = "area", bw=0.08,position = "identity",
               fill=color2use,color=color2use,alpha=0.5,adjust=1)+
  ylim(c(0,3))+
  geom_vline(xintercept=0,linetype='dashed')+
  xlab('Correlations')+
  ylab('Density')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank())
return(p)
}
#cor2plot=cor(glist_expr_df,brain_df,method='pearson')
#plot_cor(cor2plot,'#33638DFF')
  
  
# make brain plot function
plot_brain <- function(df2plot,col2plot,show.legend=F){
p=df2plot%>%mutate(label=rownames(.)) %>% 
  select(label,col2plot)%>%
  mutate(
    label=sub('.*L_','lh_',label),
    label=sub('.*R_','rh_',label),
    label=sub('_thickavg','',label))%>%
  rename('ES'=col2plot)%>%
  brain_join(dk) %>%
  reposition_brain(.~hemi+side) %>%
  ggplot(aes(fill = ES)) +
  geom_sf(show.legend = show.legend) +
  scale_fill_viridis()+
  #scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=c(-0.6,0.6),aes(title="ES"))+
  theme_void()
return(p)
}
#plot_brain(df2plot=brain_df, col2plot = 'x')


# make heatplot function
plot_heatmap<- function(df2plot){
df2heat <- df2plot%>%
  tibble::rownames_to_column() %>%
  tidyr::gather(colname, value, -rowname)
p=ggplot(df2heat, aes(y = rowname, x = colname, fill = value)) +
  geom_tile(colour="black",size=0.2)+
  #scale_x_discrete(labels = abbreviate)+
  scale_fill_viridis()+
  theme_void()+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
return(p)
}
# plot_heatmap(df2plot=brain_df)
# plot_heatmap(df2plot=glist_expr_df)


# Use Microglia as an exmaple
data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data'
gall_expr_df=get_expression.shin2stage()


celltype_df = read.csv(sprintf('%s/Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv',data_path),
                       sep = '\t',stringsAsFactors = F)
glist=celltype_df$GeneSymbol[!is.na(celltype_df$CellType) & celltype_df$CellType=='Microglia']
glist_expr_df=get_glist_expression(gall_expr_df,glist)


brain_df=read_brain_phenotype (data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data',
                               file_name='ENIGMAALC_ES.csv',
                               hem='L')

plot_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/plots/illustration'

# brain maps
p=plot_brain(brain_df,'x')
ggsave(sprintf('%s/brain_map.png',plot_path),p,height = 4, width = 6)

for (g2plot in colnames(glist_expr_df)[c(1:3,ncol(glist_expr_df))]){
  p=plot_brain(glist_expr_df,g2plot)
  ggsave(sprintf('%s/microgial_%s_brain_map.png',plot_path,g2plot),p,height = 4, width = 6)
}

####################################
# dist of correlation 
# no spin
perm.id.spin10k=readRDS(paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_spin10k.rds',sep='/')) 
perm.id.10k=perm.id.spin10k
null_brain_df=build_null_brain_df(perm.id.10k,brain_df)
rownames(null_brain_df)=rownames(brain_df)
for (brain_i in c(1:3)){

  p=plot_brain(null_brain_df,colnames(null_brain_df)[brain_i])
  ggsave(sprintf('%s/brain_map_null%s.png',plot_path,brain_i),p,height = 4, width = 6)
  
  cor_nullbrain=cor(glist_expr_df,null_brain_df[,brain_i,drop=F])
  colnames(cor_nullbrain)='x'
  p=plot_cor(cor_nullbrain,'#33638DFF')
  ggsave(sprintf('%s/microgial_cor_resample_brain_%s.png',plot_path,brain_i),p,height = 4, width = 6)
}


####################################
# dist of correlation 
# no spin
perm.id.nospin10k=readRDS(paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_nospin10k.rds',sep='/')) 
perm.id.10k=perm.id.nospin10k
null_brain_df=build_null_brain_df(perm.id.10k,brain_df)
rownames(null_brain_df)=rownames(brain_df)
for (brain_i in c(1:3)){
  
  p=plot_brain(null_brain_df,colnames(null_brain_df)[brain_i])
  ggsave(sprintf('%s/brain_map_null%s_nospin.png',plot_path,brain_i),p,height = 4, width = 6)
  
  cor_nullbrain=cor(glist_expr_df,null_brain_df[,brain_i,drop=F])
  colnames(cor_nullbrain)='x'
  p=plot_cor(cor_nullbrain,'#33638DFF')
  ggsave(sprintf('%s/microgial_cor_resample_nospin_brain_%s.png',plot_path,brain_i),p,height = 4, width = 6)
}






# heat maps
p=plot_heatmap(brain_df)
ggsave(sprintf('%s/brain_heatmap.png',plot_path),p,height = 4, width = 0.2)

p=plot_heatmap(glist_expr_df)
ggsave(sprintf('%s/microgial_heatmap.png',plot_path),p,height = 4, width = 6)


for (glist_i in c(1:3)){
  set.seed(glist_i)
  glist_resample=colnames(gall_expr_df)[sample(ncol(gall_expr_df),ncol(glist_expr_df),replace = F)]
  glist_expr_df_resample=get_glist_expression(gall_expr_df,glist_resample)
  p=plot_heatmap(glist_expr_df_resample)
  ggsave(sprintf('%s/microgial_heatmap_resample%s.png',plot_path,glist_i),p,height = 4, width = 6)
  
  cor_resample=cor(glist_expr_df_resample,brain_df,method='pearson')
  p=plot_cor(cor_resample,'#33638DFF')
  ggsave(sprintf('%s/microgial_cor_resample_gene_%s.png',plot_path,glist_i),p,height = 4, width = 6)
}


# cor
cor2plot=cor(glist_expr_df,brain_df,method='pearson')
p=plot_cor(cor2plot,'#FDE725FF')
ggsave(sprintf('%s/microgial_cor.png',plot_path),p,height = 4, width = 6)

col_grid=rgb(235, 235, 235, 100, maxColorValue = 255)

dist2plot=data.frame(x=rnorm(10000))
color2use='#33638DFF'
p1=ggplot(dist2plot, aes(x=x)) +
  stat_density(geom = "area", bw=0.1,position = "identity",
               fill=color2use,color=color2use,alpha=0.5,adjust=1)+
  geom_vline(xintercept=-1.2,color='#FDE725FF')+
  xlab('Correlations')+
  ylim(c(0,1))+
  ylab('Density')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank())+
  theme(panel.grid = element_line(color = col_grid))

dist2plot=data.frame(x=rnorm(10000))
p2=ggplot(dist2plot, aes(x=abs(x))) +
  stat_density(geom = "area", position = "identity",
               fill=color2use,color=color2use,alpha=0.5,adjust=1)+
  geom_vline(xintercept=1.2,color='#FDE725FF')+
  xlab('Correlations')+
  ylim(c(0,1))+
  xlim(c(0,3))+
  ylab('Density')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank())+
  theme(panel.grid = element_line(color = col_grid))


dist2plot=data.frame(x=rnorm(20000))
p3=ggplot(dist2plot, aes(x=abs(x)) )+
  stat_density(geom = "area", bw=0.2,position = "identity",
               fill=color2use,color=color2use,alpha=0.5,adjust=1)+
  geom_vline(xintercept=2,color='#FDE725FF')+
  xlab('Correlations')+
  ylim(c(0,2))+
  xlim(c(1,4))+
  ylab('Density')+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank())+
  theme(panel.grid = element_line(color = col_grid))

main=grid.arrange(p1,p2,p3,ncol=1)
ggsave(sprintf('%s/sigtest.png',plot_path),main,height = 8, width = 4)

