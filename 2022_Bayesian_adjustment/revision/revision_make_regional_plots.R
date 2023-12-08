## regional plots
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggseg)
library(ggridges)
library(grid)
library(gridExtra)
library(forcats)
library(HDInterval)
get_mode <- function(paramSampleVec){
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  return(mo)
}

mode_hdi <- function(x,...){
  hdi_bound=hdi(x)
  mode=get_mode(x)
  return(c(hdi_bound['lower'],mode,hdi_bound['upper']))
}
plot_study_mu <- function(all_post,
                          site_n_pad,
                          region='R_hippo',
                          xlab2show='mu',
                          ylab2show='Study (n)',
                          title2show='Hippocampus (Right)'){

df2plot=as.data.frame(all_post[[region]]) %>%
  select(-c(common_sigma,tau)) %>% 
  gather(key='Study',value='Density') %>% 
  mutate(site_n=rep(site_n_pad, each=400000)) %>% 
  arrange(site_n) %>% 
  mutate(site_idx=rep(1:22, each=400000)) %>% 
  mutate(site_y_label=paste0(site_idx,' (', site_n, ')')) %>% 
  mutate(site_y_label=as.factor(as.character(site_y_label))) 

mu_mode_hdi=mode_hdi(df2plot[df2plot$Study=='mu','Density'])


p=ggplot(df2plot) + 
  geom_density_ridges_gradient(aes(x = Density, y = site_y_label,fill = stat(quantile)),
                               quantile_lines = TRUE, quantile_fun = mode_hdi, vline_linetype = 2) +
  scale_fill_manual(values = c("skyblue4", "skyblue", "skyblue","skyblue4"), guide = "none")+
  scale_y_discrete(labels=c(as.character(unique(df2plot$site_y_label)[1:21]),'Overarching'),limits = unique(df2plot$site_y_label))+
  geom_rect(data = data.frame(x = 1),
            xmin = mu_mode_hdi[1], xmax = mu_mode_hdi[3], ymin = -Inf, ymax = Inf,
            alpha = 0.4, fill = "gray")+
  geom_vline(xintercept = mu_mode_hdi[2],linetype=2)+
  xlab(xlab2show)+
  ylab(ylab2show)+
  ggtitle(title2show)+
  theme_ridges(font_size = 11,
               font_family = "sans", # windowsFonts()-->"TT Arial"
               center_axis_labels = T)+
  theme(plot.title = element_text(hjust = 0.5))
return(p)
}


# plot raw effect size against mode value of posterior distribution
plot_raw_vs_bayes <- function(all_post,
                              all_site_es,
                              site_n,
                              region='R_hippo',
                              stat2use=c('mean','mode','median')){

  site_raw_es=all_site_es[,region]
post_df=as.data.frame(all_post[[region]]) %>%
        select(-c(common_sigma,tau,mu))

if (stat2use=='mode'){
site_bayes_es=as.numeric(apply(post_df, 2, get_mode))} 
else if(stat2use=='mean'){
  site_bayes_es=as.numeric(apply(post_df, 2, mean))}
else if(stat2use=='median'){
  site_bayes_es=as.numeric(apply(post_df, 2, median))}


mu_mode_hdi=mode_hdi(all_post[[region]][,'mu'])

df2plot=data.frame(site_bayes_es,site_raw_es) %>% 
            gather(key='stat_type',value='es') %>% 
            mutate(site_n=rep(site_n,times=2)) %>% 
            mutate(site_idx.tmp=rep(c(1:21),time=2)) %>% 
            arrange(site_n,site_idx.tmp) %>% 
            mutate(site_idx=rep(c(1:21),each=2)) %>% 
            mutate(site_y_label=sprintf('%d (%d)',site_idx,site_n))
            

p=ggplot(data=df2plot)+
  geom_point(aes(y=site_y_label,x=es,color=stat_type),alpha=0.8)+
  geom_path(data = df2plot, aes(y=site_y_label,x=es,group=site_y_label),color='black',
            arrow = arrow(length = unit(.015, "npc"), ends = "first", type = "closed" ))+
  scale_color_manual(name='',
                     labels=c('Raw', 'Bayes'),
                     values=c('firebrick1', 'deepskyblue1'))+
  theme_classic()+
  xlab('Effect size/mu')+
  ylab('Study (n)')+
  #scale_y_discrete(labels=df2plot$site_y_label,limits = factor(unique(df2plot$site_idx)))+
  scale_y_discrete(labels=as.character(unique(df2plot$site_y_label)),limits = unique(df2plot$site_y_label))+
  theme(text=element_text(family="sans",size=11), legend.position='none')+
  theme(panel.spacing = unit(1, "lines"))+
  theme(panel.grid.major.y = element_line(linetype = 'dashed'))+
  theme(panel.grid.major.x = element_line(linetype = 'solid'))+
  geom_rect(data = data.frame(x = 1),
            xmin = mu_mode_hdi[1], xmax = mu_mode_hdi[3], ymin = -Inf, ymax = Inf,
            alpha = 0.4, fill = "gray")+
  geom_vline(xintercept = mu_mode_hdi[2],linetype=2)
return(p)
}


setwd('F:/Google Drive/post-doc/Bayesian_Project/new_model/revision/no_combat')

all_site_es=read.csv('all_site_es.csv')
all_site_sd=read.csv('all_site_sd.csv')
site_n=as.numeric(readRDS('F:/Google Drive/post-doc/Bayesian_Project/new_model/site_N.rds')) #
site_n_pad=c(9999, site_n) 
# make plot for two regions
df=data.frame(region2plot=c('L_caudalmiddlefrontal_thickavg','R_lateralorbitofrontal_thickavg'),
              post2load=rep(c('Posterior_L_ct.rds','Posterior_R_ct.rds'),each=1),
              title2show=c('Caudal Middle Frontal (Left)','Lateral Orbitofrontal (Right)'))


for (region_i in c(1:2)){
all_post=readRDS(as.character(df$post2load[region_i]))
region2plot=as.character(df$region2plot[region_i])
title2show=as.character(df$title2show[region_i])

p1=plot_study_mu(all_post,
              site_n_pad,
              region=region2plot,
              xlab2show='mu',
              ylab2show='Study (n)',
              title2show=title2show)
p2=plot_raw_vs_bayes(all_post,
                      all_site_es,
                      site_n,
                      region=region2plot,
                      stat2use = 'mode')
#grid.arrange(p1,p2,ncol=1)
ggsave(sprintf('%s_study_mu.png',region2plot),p1,width = 6, height = 4,bg='white')
ggsave(sprintf('%s_study_diff.png',region2plot),p2,width = 6, height = 4,bg='white')
}

