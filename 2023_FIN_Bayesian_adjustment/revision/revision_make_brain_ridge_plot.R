
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




# Function to plot posterior distriubtions
plot_post_para <- function(all_post,
                           orig_label,
                           full_label,
                           para2plot='mu',
                           xlab2show='Effect size',
                           ylab2show='Region (Left)'){
  # key_df=read.csv('regions.csv')
  # patterns=as.character(key_df$orignal_label)
  # replacement=as.character(key_df$full_label2)
  
  # get mu column for each element of the list 
  all_post.para=lapply(all_post, function(x) x[,para2plot])
  df.para=as.data.frame(all_post.para) %>% 
    rename_all(~stringr::str_replace(.,"L_|R_","")) %>% # Rename labels
    rename_all(~stringr::str_replace_all(.,setNames(full_label,paste0(orig_label,'$')))) %>% # add $ to get exact match 
    gather(key='Region',value = 'Density') %>% 
    mutate(Regions=fct_reorder(Region,Density,.fun=get_mode,.desc = TRUE)) 
  
  
  p=ggplot(df.para) + 
    geom_density_ridges_gradient(aes(x = Density, y = Regions,fill = stat(quantile)),
                                 quantile_lines = TRUE, quantile_fun = mode_hdi, vline_linetype = 2) +
    scale_fill_manual(values = c("skyblue4", "skyblue", "skyblue","skyblue4"), guide = "none")+
    xlab(xlab2show)+
    ylab(ylab2show)+
    theme_ridges(font_size = 11,
                 font_family = "sans", # windowsFonts()-->"TT Arial"
                 center_axis_labels = T)
  return(p)
}


# Function to plot brain maps
map_post_mode <- function(all_post,
                          orig_label,
                          ggseg_label,
                          para2plot,
                          atlas2use='dk', # 'aseg'
                          hem2plot='left', # 'right'
                          title2show='Subcortical Volume (Left)',
                          limit2use=c(-0.5,0.5)){
all_post.para=lapply(all_post, function(x) x[,para2plot])
df2plot=as.data.frame(lapply(all_post.para, get_mode)) %>% 
  rename_all(~stringr::str_replace(.,"L_|R_","")) %>% # Rename labels
  #rename_all(~stringr::str_replace(.,"_thickavg|_surfavg","")) %>% # Rename labels
  rename_all(~stringr::str_replace_all(.,setNames(ggseg_label,paste0(orig_label,'$')))) %>%  # add '$' to replace the exact match
  gather(key='region',value = 'value')

p=ggseg(.data=df2plot,atlas = atlas2use,mapping=aes(fill=value), colour="darkgrey",hemisphere=hem2plot,size=.8,alpha=0.8)+
  scale_fill_gradient2(midpoint=0, low="skyblue", mid="white",high="firebrick1", space ="Lab",limits=limit2use,aes(title="M"))+
  ggtitle(title2show)+
  theme_void()+
  theme(text=element_text(size=11,  family="sans"),
        plot.title = element_text(hjust = 0.5))
return(p)
}



key_df=read.csv('F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model\\regions.csv')
orig_label=as.character(key_df$orignal_label)
full_label=as.character(key_df$full_label)
ggseg_label=as.character(key_df$ggseg_label)




setwd('F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model\\revision\\no_combat')

# CT
all_post=readRDS("Posterior_L_ct.rds")
p1=plot_post_para(all_post,orig_label,full_label, para2plot='mu',xlab2show='M', ylab2show='Regions')
p2=map_post_mode(all_post,orig_label,ggseg_label,'mu','dk','left','Cortical thickness (Left)',c(-0.3,0.3))
ggsave('Posterior_L_ct_ridge.png',p1,width = 6,height = 8,bg='white')
ggsave('Posterior_L_ct_brain.png',p2,width = 4,height = 4,bg='white')
rm(all_post,p1,p2)
gc()

all_post=readRDS("Posterior_R_ct.rds")
p1=plot_post_para(all_post,orig_label,full_label, para2plot='mu',xlab2show='M', ylab2show='Regions')
p2=map_post_mode(all_post,orig_label,ggseg_label,'mu','dk','right','Cortical thickness (Right)',c(-0.3,0.3))
ggsave('Posterior_R_ct_ridge.png',p1,width = 6,height = 8,bg='white')
ggsave('Posterior_R_ct_brain.png',p2,width = 4,height = 4,bg='white')
rm(all_post,p1,p2)
gc()
