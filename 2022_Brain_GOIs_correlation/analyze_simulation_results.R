# make plots

# plot sim results
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(dplyr)
library(nlme)




# make violin plots
prepare_df2plot_sim <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code',
                                pattern2read='all_300_sims_preserve_spatial'){
  df.nocoexp=read.csv(sprintf('%s/%s_no_coexp.csv',data_path,pattern2read))
  df.nocoexp$glist_n=as.factor(df.nocoexp$glist_n) # this is not a continuous variable
  df.nospin=read.csv(sprintf('%s/%s_ctrl_coexp_nospin.csv',data_path,pattern2read))
  df.nospin$glist_n=as.factor(df.nospin$glist_n) # this is not a continuous variable
  df.spin=read.csv(sprintf('%s/%s_ctrl_coexp_spin.csv',data_path,pattern2read))
  df.spin$glist_n=as.factor(df.spin$glist_n) # this is not a continuous variable
  
  df.nocoexp.tmp=df.nocoexp %>% 
    mutate(sim_id=paste(seed,glist_n,sep='_')) %>% # sim_id used for paired_sample test
    select(FPR.orig,FPR.abs,FPR.nsig.fdr,sim_id,glist_n,coexp) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='nocoexp') 
  
  df.nospin.tmp=df.nospin %>% 
    mutate(sim_id=paste(seed,glist_n,sep='_')) %>% 
    select(FPR.orig,FPR.abs,FPR.nsig.fdr,sim_id,glist_n,coexp) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='nospin') 
  
  df.spin.tmp=df.spin %>% 
    mutate(sim_id=paste(seed,glist_n,sep='_')) %>% 
    select(FPR.orig,FPR.abs,FPR.nsig.fdr,sim_id,glist_n,coexp) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='spin') 
  
  df2plot=do.call(rbind, list(df.nocoexp.tmp,df.nospin.tmp,df.spin.tmp))
  
  hypo_type.order=c('FPR.orig', 'FPR.abs', 'FPR.nsig.fdr')
  null_type.order=c('nocoexp','nospin','spin')
  
  df2plot$hypo_type=factor(df2plot$hypo_type,hypo_type.order)
  df2plot$null_type=factor(df2plot$null_type,null_type.order)
  return(df2plot)
}


prepare_df2plot_synGO <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code',
                                  pattern2read='all_79_sims_uniform_subset1000'){
  df.syngo=read.csv(sprintf('%s/%s.csv',data_path,pattern2read))
  df.syngo$glist_n=as.factor(df.syngo$glist_n) # this is not a continuous variable
  
  df.nospin.tmp=df.syngo %>% select(starts_with('No_spin'),glist_n,coexp) %>% 
    rename_all(~stringr::str_replace(.,"^No_spin.","")) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='nospin') 
  
  df.nocoexp.tmp=df.syngo %>% select(starts_with('No_coexp'),glist_n,coexp) %>% 
    rename_all(~stringr::str_replace(.,"^No_coexp.","")) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='nocoexp') 
  
  df.spin.tmp=df.syngo %>% select(starts_with('Spin'),glist_n,coexp) %>% 
    rename_all(~stringr::str_replace(.,"^Spin.","")) %>% 
    tidyr::gather(key='hypo_type',value='FPR',FPR.orig,FPR.abs,FPR.nsig.fdr)%>%
    mutate(null_type='spin') 
  
  df2plot=do.call(rbind, list(df.nocoexp.tmp,df.nospin.tmp,df.spin.tmp))
  hypo_type.order=c('FPR.orig', 'FPR.abs', 'FPR.nsig.fdr')
  null_type.order=c('nocoexp','nospin','spin')
  
  df2plot$hypo_type=factor(df2plot$hypo_type,hypo_type.order)
  df2plot$null_type=factor(df2plot$null_type,null_type.order)
  return(df2plot)
}


plot_sim_violin <- function(df2plot){
p=ggplot(df2plot,aes(y=FPR,fill=hypo_type,x=null_type))+
  geom_violin()+
  geom_hline(aes(yintercept=0.05),linetype="dashed")+
  xlab('Null model type')+
  ylab('Psig')+
  scale_x_discrete(labels=c('Resampling genes','Shuffling brain','Spinning brain'))+
  scale_fill_viridis(discrete = TRUE,name='Test statistic',labels=c('Average r','Average |r|','N sig'))+
  theme_minimal()+
  theme(text = element_text(size = 18))
return(p)
}


plot_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/plots'

df.sim.random=prepare_df2plot_sim(pattern2read = 'all_300_sims_totally_random' )
p=plot_sim_violin(df.sim.random)
ggsave(sprintf('%s/violin_sims_totally_random.png',plot_path),p,width = 8, height = 6)  


df.sim.sa=prepare_df2plot_sim(pattern2read ='all_300_sims_uniform_subset1000')
p=plot_sim_violin(df.sim.sa)
ggsave(sprintf('%s/violin_sims_uniform_subset1000.png',plot_path),p,width = 8, height = 6)  

df.syngo.random=prepare_df2plot_synGO(pattern2read = 'all_79_sims_totally_random')
p=plot_sim_violin(df.syngo.random)
ggsave(sprintf('%s/violin_syngo_totally_random.png',plot_path),p,width = 8, height = 6)  

df.syngo.sa=prepare_df2plot_synGO(pattern2read = 'all_79_sims_uniform_subset1000')
p=plot_sim_violin(df.syngo.sa)
ggsave(sprintf('%s/violin_syngo_uniform_subset1000.png',plot_path),p,width = 8, height = 6)  

##############################################################
### report/plot summary stats
## make bars for mean/sd/max/>0.05
get_summary_df <- function(df){
# summary_df=df %>% mutate(test_type=paste(null_type,hypo_type,sep='_')) %>% 
#   group_by(test_type) %>% 
#   summarize(min = min(FPR),
#             mean = mean(FPR),
#             max = max(FPR),
#             sd = sd(FPR),
#             perct5=sum(FPR>0.05)/n())   


summary_df=df %>% 
  group_by(null_type,hypo_type) %>% 
  summarize(min = min(FPR),
            mean = mean(FPR),
            max = max(FPR),
            sd = sd(FPR),
            perct5=sum(FPR>0.05)/n(),
            .groups = 'drop') # .groups='drop' remove group structure in df


return(summary_df)
}

plot_summary_df <- function(summary_df,
                            stat2plot='sd',
                            stat2show='Standard deviation'){
  stat2plot <- sym(stat2plot) # pass string to aes (aes_string)
 p=ggplot(summary_df,aes(x=null_type,y=!!stat2plot,fill=hypo_type))+
    geom_bar(stat="identity",position = position_dodge(0.9), size=3.5)+
    xlab('Null model type')+
    ylab(stat2show)+
    scale_x_discrete(labels=c('Resampling genes','Shuffling brain','Spinning brain'))+
    scale_fill_viridis(discrete = TRUE,name='Test statistic',labels=c('Average r','Average |r|','N sig'))+
    theme_minimal()+
    theme(legend.position="none",
          #axis.text.x=element_text(angle=90, hjust=1),
          axis.title.x=element_blank(),
          text = element_text(size = 12))
 return(p)
}


make_all_summary_plot <- function(summary_df){
  stat2plot_set=c('sd','max')
  stat2show_set=c('Standard deviation','Maximum')
  p_list=list()
  for (si in c(1:length(stat2plot_set))){
    p=plot_summary_df (summary_df,
                       stat2plot_set[si],
                       stat2show_set[si])
    p_list[[si]]=p
  }
  
  p_main=grid.arrange(grobs=p_list,ncol=length(stat2plot_set))
  return(p_main)
}

table_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/tables'


df.sim.random.sum=get_summary_df(df.sim.random)
p=make_all_summary_plot(df.sim.random.sum)
ggsave(sprintf('%s/summary_bar_sims_totally_random.png',plot_path),p,width = 8, height = 4)
write.csv(df.sim.random.sum,sprintf('%s/summary_bar_sims_totally_random.csv',table_path),row.names = F)

df.sim.sa.sum=get_summary_df(df.sim.sa)
p=make_all_summary_plot(df.sim.sa.sum)
ggsave(sprintf('%s/summary_bar_sims_uniform_subset1000.png',plot_path),p,width = 8, height = 4)  
write.csv(df.sim.sa.sum,sprintf('%s/summary_bar_sims_uniform_subset1000.csv',table_path),row.names = F)

df.syngo.random.sum=get_summary_df(df.syngo.random)
p=make_all_summary_plot(df.syngo.random.sum)
ggsave(sprintf('%s/summary_bar_syngo_totally_random.png',plot_path),p,width = 8, height = 4)
write.csv(df.syngo.random.sum,sprintf('%s/summary_bar_syngo_totally_random.csv',table_path),row.names = F)

df.syngo.sa.sum=get_summary_df(df.syngo.sa)
p=make_all_summary_plot(df.syngo.sa.sum)
ggsave(sprintf('%s/summary_bar_syngo_uniform_subset1000.png',plot_path),p,width = 8, height = 4)
write.csv(df.syngo.sa.sum,sprintf('%s/summary_bar_syngo_uniform_subset1000.csv',table_path),row.names = F)

####
####do F test in variance and mean and make a heatmap
#df=df.sim.sa
caculate_pairwise_F <- function(df){
  df=df %>% mutate(test_type=paste(null_type,hypo_type,sep='_')) 
  
  test_type.unique=unique(df$test_type)
  test_vals =matrix(nrow = length(test_type.unique), ncol = length(test_type.unique))
  p_vals =matrix(nrow = length(test_type.unique), ncol = length(test_type.unique))
  for (i in (1:length(test_type.unique))) {
    for (j in (i:length(test_type.unique))){
      if (i==j){test_vals[i,j]=NA
        p_vals[i,j]=NA}
      else{
       
        #Fligner-Killeen test may be more appropriate
        x=df %>% filter(test_type==test_type.unique[i]) %>% select(FPR) %>% mutate(group=2)
        y=df %>% filter(test_type==test_type.unique[j]) %>% select(FPR) %>% mutate(group=1)
        xy=rbind(x,y)
        xy$group=factor(xy$group,levels = c(1,2))
        tmp=var.test(FPR ~ group, data = xy) 
        #tmp=fligner.test(FPR ~ group, data = xy)
        test_vals[i,j] = tmp$statistic
        p_vals[i,j] = tmp$p.value }
      }
  }
    
    test_vals[lower.tri(test_vals)]=t(test_vals)[lower.tri(test_vals)]
    p_vals[lower.tri(p_vals)]=t(p_vals)[lower.tri(p_vals)]

    labels2show=c('Resampling genes-Mean r',
                  'Resampling genes-Mean |r|',
                  'Resampling genes-N sig',
                  'Shuffling brain-Mean r',
                  'Shuffling brain-Mean |r|',
                  'Shuffling brain-N sig',
                  'Spinning brain-Mean r',
                  'Spinning brain-Mean |r|',
                  'Spinning brain-N sig')
    dimnames(test_vals)=list(labels2show,labels2show)
    dimnames(p_vals) =list(labels2show,labels2show)
    result=list()
    result[['test_vals']]=test_vals
    result[['p_vals']]=p_vals
    return(result)
  }

library(corrplot)
cor_plot <- function(res,...){
  # type: upper or lower 
  col2use=colorRampPalette(c("firebrick1","steelblue1"))
  corrplot(res$test_vals, method="color", col=col2use(50),  
           addCoef.col = "black", 
           tl.col="black", tl.srt=45, 
           p.mat = res$p_vals, sig.level = 0.05/36, insig = "blank", #Bonf correction
           diag=FALSE ,
           number.cex=0.8,cl.ratio = 0.3, addgrid.col='grey',...
  )
}
  
sim.random.res=caculate_pairwise_F(df.sim.random)
png(file=sprintf('%s/sim_totally_random_Ftest.png',plot_path),
    width=1200, height=1200,res = 150)
cor_plot(sim.random.res,is.corr =F,type="upper")
dev.off()

sim.sa.res=caculate_pairwise_F(df.sim.sa)

png(file=sprintf('%s/sim_uniform_subset1000_Ftest.png',plot_path),
    width=1200, height=1200,res = 150)
cor_plot(sim.sa.res,is.corr =F,type="upper")
dev.off()

syngo.random.res=caculate_pairwise_F(df.syngo.random)
png(file=sprintf('%s/syngo_totally_random_Ftest.png',plot_path),
    width=1200, height=1200,res = 150)
cor_plot(syngo.random.res,is.corr =F,type="upper")
dev.off()

syngo.sa.res=caculate_pairwise_F(df.syngo.sa)
png(file=sprintf('%s/syngo_uniform_subset_Ftest.png',plot_path),
    width=1200, height=1200,res = 150)
cor_plot(syngo.sa.res,is.corr =F,type="upper")
dev.off()










################################################################################
#####plot scatter plots


plot_sim_results <- function(df2plot, col2plot='FPR.orig',ylim2show=c(-0.05,0.7)){
  p=ggplot(df2plot,aes(coexp, df2plot[,col2plot], colour=glist_n))+
    geom_point(position = position_jitter(w = 0.01, h = 0.01),
               alpha=0.6,
               size=2)+
    scale_color_viridis(discrete = TRUE)+
    geom_hline(yintercept=0.05,linetype="dashed", color = "darkgrey",size=0.8) +
    #annotate(geom='text',x=-0.2,y=0.05,label = 'FPR=5%', vjust = -2.4,color='black',angle=45,size=3.2)+
    ylim(ylim2show)+
    #xlab('Co-expression')+
    xlim(c(-0.25,0.8))+
    #ylab('False Positive Rate')+
    #ggtitle("Test r values")+
    theme_minimal()+
    theme(legend.position = "none",axis.title.x=element_blank(),axis.title.y=element_blank())
  return(p)  
}

### Figure S2-Type1-Maps
# plot for totally_random
df.nocoexp=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_totally_random_no_coexp.csv')
df.nocoexp$glist_n=as.factor(df.nocoexp$glist_n) 
p1=plot_sim_results(df.nocoexp,'FPR.orig',c(-0.05,0.7))
p2= plot_sim_results(df.nocoexp,'FPR.abs',c(-0.05,0.2))
p3=plot_sim_results(df.nocoexp,'FPR.nsig',c(-0.05,0.2))
p4=plot_sim_results(df.nocoexp,'FPR.nsig.fdr',c(-0.05,0.2))


df.nospin=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_totally_random_ctrl_coexp_nospin.csv')
df.nospin$glist_n=as.factor(df.nospin$glist_n) 
p5=plot_sim_results(df.nospin,'FPR.orig',c(-0.05,0.15))
p6=plot_sim_results(df.nospin,'FPR.abs',c(-0.05,0.15))
p7=plot_sim_results(df.nospin,'FPR.nsig',c(-0.05,0.15))
p8=plot_sim_results(df.nospin,'FPR.nsig.fdr',c(-0.05,0.15))


df.spin=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_totally_random_ctrl_coexp_spin.csv')
df.spin$glist_n=as.factor(df.spin$glist_n) 
p9=plot_sim_results(df.spin,'FPR.orig',c(-0.05,0.15))
p10=plot_sim_results(df.spin,'FPR.abs',c(-0.05,0.15))
p11=plot_sim_results(df.spin,'FPR.nsig',c(-0.05,0.15))
p12=plot_sim_results(df.spin,'FPR.nsig.fdr',c(-0.05,0.15))

# main=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4)
main=grid.arrange(p1,p2,p4,p5,p6,p8,p9,p10,p12,ncol=3)
ggsave(sprintf('%s/scatter_coexp_sim_results_totally_random.png',plot_path),main,width = 12, height = 8)


### Figure S2-Type2-Maps
# uniform_subset1000
df.nocoexp=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_uniform_subset1000_no_coexp.csv')
df.nocoexp$glist_n=as.factor(df.nocoexp$glist_n) 
p1=plot_sim_results(df.nocoexp,'FPR.orig',c(-0.05,0.7))
p2= plot_sim_results(df.nocoexp,'FPR.abs',c(-0.05,0.6))
p3=plot_sim_results(df.nocoexp,'FPR.nsig',c(-0.05,0.3))
p4=plot_sim_results(df.nocoexp,'FPR.nsig.fdr',c(-0.05,0.3))

df.nospin=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_uniform_subset1000_ctrl_coexp_nospin.csv')
df.nospin$glist_n=as.factor(df.nospin$glist_n) 
p5=plot_sim_results(df.nospin,'FPR.orig',c(-0.05,0.6))
p6=plot_sim_results(df.nospin,'FPR.abs',c(-0.05,0.6))
p7=plot_sim_results(df.nospin,'FPR.nsig',c(-0.05,0.6))
p8=plot_sim_results(df.nospin,'FPR.nsig.fdr',c(-0.05,0.6))

df.spin=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_300_sims_uniform_subset1000_ctrl_coexp_spin.csv')
df.spin$glist_n=as.factor(df.spin$glist_n) 
p9=plot_sim_results(df.spin,'FPR.orig',c(-0.05,0.3))
p10=plot_sim_results(df.spin,'FPR.abs',c(-0.05,0.3))
p11=plot_sim_results(df.spin,'FPR.nsig',c(-0.05,0.3))
p12=plot_sim_results(df.spin,'FPR.nsig.fdr',c(-0.05,0.3))

# main=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4)
main=grid.arrange(p1,p2,p4,p5,p6,p8,p9,p10,p12,ncol=3)
ggsave(sprintf('%s/scatter_coexp_sim_results__uniform_subset1000.png',plot_path),main,width = 12, height = 8)





### Figure S4-Type1-Maps
# syn GO + totally random brain map
df.syngo=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_79_sims_totally_random.csv')
df.syngo$glist_n=as.factor(df.syngo$glist_n)
p1=plot_sim_results(df.syngo,'No_coexp.FPR.orig',c(-0.05,0.7))
p2= plot_sim_results(df.syngo,'No_coexp.FPR.abs',c(-0.05,0.3))
p3=plot_sim_results(df.syngo,'No_coexp.FPR.nsig',c(-0.05,0.3))
p4=plot_sim_results(df.syngo,'No_coexp.FPR.nsig.fdr',c(-0.05,0.3))
p5=plot_sim_results(df.syngo,'No_spin.FPR.orig',c(-0.05,0.3))
p6= plot_sim_results(df.syngo,'No_spin.FPR.abs',c(-0.05,0.3))
p7=plot_sim_results(df.syngo,'No_spin.FPR.nsig',c(-0.05,0.3))
p8=plot_sim_results(df.syngo,'No_spin.FPR.nsig.fdr',c(-0.05,0.3))

p9=plot_sim_results(df.syngo,'Spin.FPR.orig',c(-0.05,0.3))
p10= plot_sim_results(df.syngo,'Spin.FPR.abs',c(-0.05,0.3))
p11=plot_sim_results(df.syngo,'Spin.FPR.nsig',c(-0.05,0.3))
p12=plot_sim_results(df.syngo,'Spin.FPR.nsig.fdr',c(-0.05,0.3))


# main=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4)
main=grid.arrange(p1,p2,p4,p5,p6,p8,p9,p10,p12,ncol=3)
ggsave(sprintf('%s/scatter_coexp_synGO_results_totally_random.png',plot_path),main,width = 12, height = 8)

### Figure S4-Type2-Maps
# syn GO + uniform_subset1000
df.syngo=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/all_79_sims_uniform_subset1000.csv')
df.syngo$glist_n=as.factor(df.syngo$glist_n)
p1=plot_sim_results(df.syngo,'No_coexp.FPR.orig',c(-0.05,0.7))
p2= plot_sim_results(df.syngo,'No_coexp.FPR.abs',c(-0.05,0.5))
p3=plot_sim_results(df.syngo,'No_coexp.FPR.nsig',c(-0.05,0.5))
p4=plot_sim_results(df.syngo,'No_coexp.FPR.nsig.fdr',c(-0.05,0.5))
p5=plot_sim_results(df.syngo,'No_spin.FPR.orig',c(-0.05,0.5))
p6= plot_sim_results(df.syngo,'No_spin.FPR.abs',c(-0.05,0.5))
p7=plot_sim_results(df.syngo,'No_spin.FPR.nsig',c(-0.05,0.5))
p8=plot_sim_results(df.syngo,'No_spin.FPR.nsig.fdr',c(-0.05,0.6))

p9=plot_sim_results(df.syngo,'Spin.FPR.orig',c(-0.05,0.3))
p10= plot_sim_results(df.syngo,'Spin.FPR.abs',c(-0.05,0.3))
p11=plot_sim_results(df.syngo,'Spin.FPR.nsig',c(-0.05,0.3))
p12=plot_sim_results(df.syngo,'Spin.FPR.nsig.fdr',c(-0.05,0.3))


# main=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4)
main=grid.arrange(p1,p2,p4,p5,p6,p8,p9,p10,p12,ncol=3)
ggsave(sprintf('%s/scatter_coexp_synGO_results_all_79_sims_uniform_subset1000.png',plot_path),main,width = 12, height = 8)



# do correlation between co-expression and Psig
library(lmerTest)
library(purrr)
library(broom)
library(broom.mixed)
library(tidyr)
# # LME model with random slope and intercept 
run_lme <- function(df){
  model=lmer(FPR~coexp+(0+coexp|glist_n)+(1|glist_n),df)
  return(model)
}

run_lm <- function(df){
  model=lm(FPR~coexp,df)
}


get_stat_table <- function(df){
  df.res=df %>%
    group_by(null_type,hypo_type) %>% 
    nest() %>% 
    mutate(l_models=map(data,run_lm)) %>%  # run_lm or run_lme
    mutate(tidy_db =  map(l_models,tidy)) %>% # lme_models %>% map(glance)
    unnest(tidy_db) %>% 
    filter(term=='coexp') %>% 
    select(hypo_type,null_type,statistic,p.value)
  return(df.res)
}

### Table S2 
table_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/tables'
df.sim.random=prepare_df2plot_sim(pattern2read = 'all_300_sims_totally_random' )
df.sim.random.res = get_stat_table(df.sim.random) 

df.sim.sa=prepare_df2plot_sim(pattern2read ='all_300_sims_uniform_subset1000')
df.sim.sa.res = get_stat_table(df.sim.sa)
write.csv(df.sim.random.res,sprintf('%s/cor_coexp_sims_totally_random.csv',table_path),row.names = F)
write.csv(df.sim.sa.res,sprintf('%s/cor_coexp_sims_uniform_subset1000.csv',table_path),row.names = F)

### Table S6
df.syngo.random=prepare_df2plot_synGO(pattern2read = 'all_79_sims_totally_random')
df.syngo.random.res = get_stat_table(df.syngo.random) 
  
df.syngo.sa=prepare_df2plot_synGO(pattern2read = 'all_79_sims_uniform_subset1000')
df.syngo.sa.res=get_stat_table(df.syngo.sa)
write.csv(df.syngo.random.res,sprintf('%s/cor_coexp_syngo_totally_random.csv',table_path),row.names = F)
write.csv(df.syngo.sa.res,sprintf('%s/cor_coexp_syngo_uniform_subset1000.csv',table_path),row.names = F)

# plot fitted lines
 # ggplot(fortify(model), aes(coexp, FPR.orig, color=glist_n)) +
 #  geom_point(position = position_jitter(w = 0.01, h = 0.01),
 #    alpha=0.4, size=2)+scale_color_viridis(discrete = TRUE)+
 #  stat_summary(aes(y=.fitted), fun=mean, geom="line",linetype="dashed",size=1)+# connect fitted y values by line
 #  xlab('Co-expression')+
 #  ylab('P(Sig)')+
 #  theme_minimal()+
 #  theme(legend.position = "none",
 #        axis.title.x=element_blank(),
 #        axis.title.y=element_blank())

# 
#  model=lmer(FPR.abs~coexp+(0+coexp|glist_n),df.nocoexp) 
#  summary(model)
#  ggplot(fortify(model), aes(coexp, FPR.abs, color=glist_n)) +
#    geom_point(position = position_jitter(w = 0.01, h = 0.01),
#               alpha=0.4, size=2)+scale_color_viridis(discrete = TRUE)+
#    stat_summary(aes(y=.fitted), fun=mean, geom="line",linetype="dashed",size=1)+# connect fitted y values by line
#    #xlab('Co-expression')+
#    #ylab('False Positive Rate')+
#    theme_minimal()+
#    theme(legend.position = "none",
#          axis.title.x=element_blank(),
#          axis.title.y=element_blank())
#  
#  model=lmer(FPR.nsig.fdr~coexp+(0+coexp|glist_n),df.nocoexp) 
#  summary(model)
#  ggplot(fortify(model), aes(coexp, FPR.nsig.fdr, color=glist_n)) +
#    geom_point(position = position_jitter(w = 0.01, h = 0.01),
#               alpha=0.4, size=2)+scale_color_viridis(discrete = TRUE)+
#    stat_summary(aes(y=.fitted), fun=mean, geom="line",linetype="dashed",size=1)+# connect fitted y values by line
#    #xlab('Co-expression')+
#    #ylab('False Positive Rate')+
#    theme_minimal()+
#    theme(legend.position = "none",
#          axis.title.x=element_blank(),
#          axis.title.y=element_blank())
 

  

# BayesFactor gave similar results. So report t
# library(BayesFactor)
# 
# summary(lmer(FPR.orig~coexp+(1|glist_n),df.nocoexp))
# 
# full_BF = lmBF(FPR.orig ~ coexp + glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.orig~ glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.orig ~ coexp + glist_n, data = df.nospin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.orig~ glist_n, data = df.nospin, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.orig ~ coexp + glist_n, data = df.spin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.orig~ glist_n, data = df.spin, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# # abs
# full_BF = lmBF(FPR.abs ~ coexp + glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.abs~ glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.abs ~ coexp + glist_n, data = df.nospin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.abs~ glist_n, data = df.nospin, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.abs ~ coexp + glist_n, data = df.spin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.abs~ glist_n, data = df.spin, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# 
# # fdr.n
# full_BF = lmBF(FPR.nsig.fdr ~ coexp + glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.nsig.fdr~ glist_n, data = df.nocoexp, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.nsig.fdr ~ coexp + glist_n, data = df.nospin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.nsig.fdr~ glist_n, data = df.nospin, whichRandom = 'glist_n')
# full_BF / null_BF
# 
# full_BF = lmBF(FPR.nsig.fdr ~ coexp + glist_n, data = df.spin, whichRandom = 'glist_n')
# null_BF = lmBF(FPR.nsig.fdr~ glist_n, data = df.spin, whichRandom = 'glist_n')
# full_BF / null_BF



  
  
  

  
  
  
  
  
  
  
