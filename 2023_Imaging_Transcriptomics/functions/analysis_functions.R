library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(lmerTest)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(pheatmap)
library(ggplotify)

# get pvals information from res.df
# nest_by==geneSet: pvals.nested is grouped by geneSet. then psig is the prob of detecting sig geneSet in geneSetList for each brain
# nest_by==brain:  pvals.nested is grouped by brain, then psig is the prob of detecting sig brain in all brain maps for each geneSet
# heat_plot: generate heat_plot for the significant association (not recommended for large geneSetList)
# caculate_sim: caculate the similarity for the significant association across test statstics
# sim.method: what measure will be used to measure the similarity aross test statstics, default is jaccard
get_pvals_nested <- function(res.df,
                             threshold = 0.05,
                             nest_by=c('geneSet','brain'),
                             sim.method=c("jaccard","matching","pearson","phi","gupta",
                                          "affinity","cosine","euclidean","dice"),
                             heat_plot=F,
                             caculate_sim=T){
  
  
  nest_by=match.arg(nest_by)
  sim.method=match.arg(sim.method)
  
  
  # select start with pvals and create adj.pvals
  pvals.df=res.df %>%
    select(c(starts_with('pvals'), null_brain,geneSet)) %>% 
    mutate(brain_id=paste0('brain_',null_brain)) %>% 
    select(-null_brain)
  
  if (nest_by=='brain'){
    pvals.nested = pvals.df %>% #need to keep geneSet for the following
      group_by(brain_id) %>%
      nest() %>% 
      mutate(pvals=map(data, function(x) {tibble::column_to_rownames(x,'geneSet')})) %>% 
      select(-data)%>% 
      mutate(pvals.fdr = map(pvals, function(x) {apply(x, 2, p.adjust, method = 'fdr')}))
  } else if (nest_by=='geneSet'){
    pvals.nested = pvals.df %>% 
      group_by(geneSet) %>%
      nest() %>% 
      mutate(pvals=map(data, function(x) {x %>% tibble::column_to_rownames('brain_id')})) %>% 
      select(-data) %>% 
      mutate(pvals.fdr = map(pvals, function(x) {apply(x, 2, p.adjust, method = 'fdr')}))
  }
  if (heat_plot){
    pvals.nested=pvals.nested %>% 
      mutate(pvals.plot=map(pvals, heat_plot_pvals_df,threshold=threshold)) %>% 
      mutate(pvals.fdr.plot=map(pvals.fdr,heat_plot_pvals_df,threshold=threshold))
  }
  if (caculate_sim){
    pvals.nested=pvals.nested %>% 
      mutate(pvals.sim=map(pvals, get_pvals_df_similarity,threshold=threshold,method=sim.method)) %>% 
      mutate(pvals.fdr.sim=map(pvals.fdr, get_pvals_df_similarity,threshold=threshold,method=sim.method))
  }
  return(pvals.nested)
}

# heat plot for the pvals.df
heat_plot_pvals_df <- function(pvals.df,
                               threshold=0.05){
  sig.df=pvals.df<threshold
  sig.df=apply(sig.df,2,as.numeric)
  rownames(sig.df)=rownames(pvals.df)
  p=reshape2::melt(sig.df) %>% 
    rename(Statistic=Var2) %>% 
    mutate(stat_type=stringr::str_replace(Statistic,'pvals.','')) %>% 
    mutate(stat_sort=factor(stat_type, level=c('mean',
                                               'meanabs',
                                               'meansqr',
                                               'maxmean',
                                               'median',
                                               'sign_test',
                                               'sig_n',
                                               'rank_sum',
                                               'ks_orig',
                                               'ks_weighted'))) %>% 
    ggplot(aes(x=stat_sort,y=Var1,fill=as.factor(value)))+
    geom_tile(color = "white",
              lwd = 0.2,
              linetype = 1) +
    coord_fixed()+
    scale_fill_manual(values = c('steelblue1','firebrick1'))+
    scale_y_discrete(limits=rev)+
    ylab('GeneSet or Brain')+
    scale_x_discrete(name='Test statistic',
                     labels=c('Mean',
                              'Meanabs',
                              'Meansqr',
                              'Maxmean',
                              'Median',
                              'Sign Test',
                              'Sig Number',
                              'Wilcoxon Rank Sum',
                              'KS',
                              'Weighted KS'),
                     position = 'top')+
    theme_minimal()+
    theme(axis.text.x = element_text(size=12,angle=45,hjust=0.05,vjust=0.5),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14,face="bold"),
          legend.position="none")
  return(p)
}


heat_plot2_pvals_df <- function(pvals.df,
                                main2show,
                                threshold=0.05){
  sig.df=pvals.df<threshold
  sig.df=apply(sig.df,2,as.numeric)
  rownames(sig.df)=rownames(pvals.df)
  colnames(sig.df)=recode(stringr::str_remove(colnames(sig.df),'pvals.'),
                          mean='Mean',
                          median='Median',
                          meanabs='Meanabs',
                          meansqr='Meansqr',
                          maxmean='Maxmean',
                          sig_n='Sig Number',
                          #sign_test='Sign Test', 
                          #rank_sum='Wilcoxon Rank Sum',
                          ks_orig='Kolmogorov-Smirnov (KS)',
                          ks_weighted='Weighted KS')
  sig.df=sig.df[,-which(colnames(sig.df) %in% c("sign_test", "rank_sum"))]
  
  pheat=pheatmap(sig.df,
                 color = colorRampPalette(c('#8ebcd9','#fdcb82'))(100),
                 cellwidth = 20,
                 cellheight = 15,
                 angle_col = 45,
                 fontsize = 14,
                 na_col='grey30',
                 cluster_cols = F,
                 cluster_rows = F,
                 legend = F,
                 main=main2show,
                 show_rownames = FALSE)
  return(pheat)
}



get_pvals_df_similarity <- function(pvals.df, 
                                    threshold=0.05,
                                    method=c("matching","jaccard","pearson","phi","gupta",
                                             "affinity","cosine","euclidean","dice")){
  method=match.arg(method)
  sig.df=pvals.df<threshold
  mat=apply(sig.df,2,as.numeric)
  rownames(sig.df)=rownames(pvals.df)
  present = !is.na(mat) & mat > 0
  sim=1-as.matrix(arules::dissimilarity(t(present), method = method))
  return(sim)
}

# get psig for the pvals.nested structure
# if pvals.nested is grouped by brain. then psig is the prob of detecting sig geneSet in geneSetList for each brain
# if pvals.nested is grouped by geneSet, then psig is the prob of detecting sig brain in all brain maps for each geneSet
# if_fdr=T: return the results based on fdr corrected p values. 
get_psig <- function(pvals.nested,
                     if_fdr=F){ 
  if (!if_fdr){
  psig = pvals.nested %>%
    mutate(psig = map(pvals, function(x) {
      colMeans(x < 0.05)
    })) %>%
    select(psig) %>%
    unnest_wider(psig) %>% # use unnest_wider
    ungroup() %>%
    rename_all( ~ stringr::str_replace(., '^pvals.', ''))
  } else {
  psig = pvals.nested %>%
    mutate(psig = map(pvals.fdr, function(x) {colMeans(x < 0.05)})) %>%
    select(psig) %>%
    unnest_wider(psig) %>% # use unnest_wider
    ungroup() %>%
    rename_all( ~ stringr::str_replace(., '^pvals.', ''))
  }
  return(psig)
}

# get mean and std/se for the psig
summarise_psig <- function(psig){
  s_psig=psig %>% 
    pivot_longer(-1, names_to = "stat", values_to = "value") %>% 
    group_by(stat) %>% 
    summarize(mean_val=mean(value),
              sd_val=sd(value),
              n_val=n(),
              se_val=sd_val/sqrt(n_val),
              upper_limit=mean_val+se_val,
              lower_limit=mean_val-se_val
    )
  return(s_psig)
}

# get the avearge simarity matrix for all null_brain/geneSets
# if_fdr=T: return the results based on fdr corrected p values. 
get_mean_sim <- function(pvals.nested,
                         if_fdr=F){
  if (!if_fdr){
  sim.list=pvals.nested$pvals.sim
  } else (
  sim.list=pvals.nested$pvals.fdr.sim
  )
  sim.avg=Reduce("+", sim.list) / length(sim.list)
  return(sim.avg)
}

heatmap_sim.avg <- function(sim.avg,
                            cutree=2,
                            title2show='test'){
  sim.avg=sim.avg[!rownames(sim.avg) %in% c('pvals.sign_test','pvals.rank_sum'), ] 
  sim.avg=sim.avg[,!colnames(sim.avg) %in% c('pvals.sign_test','pvals.rank_sum')] 
  # recode sim.mx col/row names
  colnames(sim.avg)=recode(stringr::str_remove(colnames(sim.avg),'pvals.'),
                          mean='Mean',
                          median='Median',
                          meanabs='Meanabs',
                          meansqr='Meansqr',
                          maxmean='Maxmean',
                          sig_n='Sig Number',
                          #sign_test='Sign Test', 
                          #rank_sum='Wilcoxon Rank Sum',
                          ks_orig='KS',
                          ks_weighted='Weighted KS')
  
  rownames(sim.avg)=recode(stringr::str_remove(rownames(sim.avg),'pvals.'),
                          mean='Mean',
                          median='Median',
                          meanabs='Meanabs',
                          meansqr='Meansqr',
                          maxmean='Maxmean',
                          sig_n='Sig Number',
                          #sign_test='Sign Test', 
                          #rank_sum='Wilcoxon Rank Sum',
                          ks_orig='KS',
                          ks_weighted='Weighted KS')
  
  diag(sim.avg)=NA
  p=pheatmap(sim.avg,
           angle_col = 45,
           cutree_rows = cutree, 
           cutree_cols =cutree,
           fontsize = 14,
           display_numbers = T,
           number_format = "%.2f",
           na_col='grey30',
           main=title2show)
  
  p$gtable$grobs[[1]]
  # The result of the operation is £ºtext[GRID.text. Numbers ] 
  p$gtable$grobs[[1]]$label
  p$gtable$grobs[[1]]$x
  # Example heat map code £º
  p$gtable$grobs[[1]]$gp$fontsize <- 16
  p$gtable$grobs[[1]]$x <- unit(-0.2,'npc')
  p$gtable$grobs[[1]]$y <- unit(0.3,'npc')
  p$gtable$grobs[[1]]$just <- c('left','bottom')
  grid.newpage()
  p
  return(p)
}


save_pheatmap_png <- function(x, 
                              filename, 
                              width=1200, 
                              height=1000, 
                              res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}




# plot violin plot for the raw scatters
# psig.list is the list of psig for each null brain/geneSet
# show_shuffle: if show the results of shuffling brain
plot_raw_psig_list <- function(psig.list,
                               show_shuffle=F,
                               ylab2show='Psig',# results of shuffling brain was shown as SI
                               title2show='Psig',
                               title_adj=-0.1){# 0.5 is to center the title
  
  # plot mutiple psig in the list
  psig.list=lapply(seq_along(psig.list), function(i){psig.list[[i]] %>% mutate(null_type=names(psig.list)[[i]])})
  psig2plot=do.call(rbind,psig.list)
  
  
  if (show_shuffle){
    p =psig2plot %>% 
      select(-c(1)) %>% 
      pivot_longer(-null_type,names_to = "stat", values_to = "value") %>%
      mutate(null_type_sort=factor(null_type, level = c('random_gene',
                                                        'random_brain',
                                                        'spin_brain'))) %>% 
      filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
      mutate(stat_sort=factor(stat, level=c('mean',
                                            'meanabs',
                                            'meansqr',
                                            'maxmean',
                                            'median',
                                            #'sign_test',
                                            'sig_n',
                                            #'rank_sum',
                                            'ks_orig',
                                            'ks_weighted'))) %>% 
      mutate(vlineposition=as.numeric(factor(null_type, levels = c('random_gene',
                                                                   'random_brain',
                                                                   'spin_brain'))),
             vlineposition=vlineposition+.5,
             vlineposition=ifelse(vlineposition == max(vlineposition), NA, vlineposition)) %>% 
      ggplot(aes(x=null_type_sort,y=value,fill=stat_sort))+
      geom_violin()+
      ylab(ylab2show)+
      xlab('Null model type')+
      geom_hline(aes(yintercept=0.05),linetype="dashed")+
      ggtitle(title2show)+
      geom_vline(aes(xintercept = vlineposition),color='grey',alpha=0.5)+
      scale_x_discrete(labels= c('Resampling gene',
                                 'Shuffling brain',
                                 'Spinning brain'))+
      scale_fill_viridis(discrete = TRUE,
                         name='Test statistic',
                         labels=c('Mean',
                                  'Meanabs',
                                  'Meansqr',
                                  'Maxmean',
                                  'Median',
                                  #'Sign Test',
                                  'Sig Number',
                                  #'Wilcoxon Rank Sum',
                                  'KS',
                                  'Weighted KS'))+
      theme_minimal(base_size = 20)+
      theme(axis.text.x = element_text(vjust=10),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust =title_adj)) 
    
  } else {
    p =psig2plot %>% 
      select(-c(1)) %>% 
      pivot_longer(-null_type,names_to = "stat", values_to = "value") %>%
      filter(!null_type=='random_brain') %>% 
      mutate(null_type_sort=factor(null_type, level = c('random_gene',
                                                        'spin_brain'))) %>%
      filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
      mutate(stat_sort=factor(stat, level=c('mean',
                                            'meanabs',
                                            'meansqr',
                                            'maxmean',
                                            'median',
                                            #'sign_test',
                                            'sig_n',
                                            #'rank_sum',
                                            'ks_orig',
                                            'ks_weighted'))) %>% 
      mutate(vlineposition=as.numeric(factor(null_type, levels = c('random_gene',
                                                                   'spin_brain'))),
             vlineposition=vlineposition+.5,
             vlineposition=ifelse(vlineposition == max(vlineposition), NA, vlineposition)) %>% 
      ggplot(aes(x=null_type_sort,y=value,fill=stat_sort))+
      geom_violin()+
      ylab(ylab2show)+
      xlab('Null model type')+
      geom_hline(aes(yintercept=0.05),linetype="dashed")+
      ggtitle(title2show)+
      geom_vline(aes(xintercept = vlineposition),color='grey',alpha=0.5)+
      scale_x_discrete(labels= c('Competetive null models',
                                 'Self-contained null models'))+
      scale_fill_viridis(discrete = TRUE,
                         name='Test statistic',
                         labels=c('Mean',
                                  'Meanabs',
                                  'Meansqr',
                                  'Maxmean',
                                  'Median',
                                  #'Sign Test',
                                  'Sig Number',
                                  #'Wilcoxon Rank Sum',
                                  'KS',
                                  'Weighted KS'))+
      theme_minimal(base_size = 20)+
      theme(axis.text.x = element_text(vjust=10),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust =title_adj))   
  }
  
  return(p)
}


# plot summaries (bar+se) of the scatters
plot_summary_psig_list <- function(psig.list,
                                   show_shuffle=F, 
                                   ylab2show='Psig',
                                   title2show='Psig',
                                   title_adj=-0.1){
  summary.list=lapply(psig.list, summarise_psig)
  summary.list=lapply(seq_along(summary.list), function(i){summary.list[[i]] %>% mutate(null_type=names(summary.list)[[i]])})
  summary2plot=do.call(rbind,summary.list) 
  
  if (show_shuffle){
  p =summary2plot  %>% 
    mutate(null_type_sort=factor(null_type, level = c('random_gene',
                                                      'random_brain',
                                                      'spin_brain'))) %>% 
    filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
    mutate(stat_sort=factor(stat, level=c('mean',
                                          'meanabs',
                                          'meansqr',
                                          'maxmean',
                                          'median',
                                          #'sign_test',
                                          'sig_n',
                                          #'rank_sum',
                                          'ks_orig',
                                          'ks_weighted'))) %>% 
    mutate(vlineposition=as.numeric(factor(null_type, levels = c('random_gene',
                                                                 'random_brain',
                                                                 'spin_brain'))),
           vlineposition=vlineposition+.5,
           vlineposition=ifelse(vlineposition == max(vlineposition), NA, vlineposition)) %>% 
    ggplot(aes(x=null_type_sort,y=mean_val,fill=stat_sort))+
    geom_bar(stat="identity",position="dodge2") +
    geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit),position="dodge2")+
    ylab(ylab2show)+
    xlab('Null model type')+
    geom_hline(aes(yintercept=0.05),linetype="dashed")+
    ggtitle(title2show)+
    geom_vline(aes(xintercept = vlineposition),color='grey',alpha=0.5)+
    scale_x_discrete(labels= c('Resampling gene',
                               'Shuffling brain',
                               'Spinning brain'))+
    scale_fill_viridis(discrete = TRUE,
                       name='Test statistic',
                       labels=c('Mean',
                                'Meanabs',
                                'Meansqr',
                                'Maxmean',
                                'Median',
                                #'Sign Test',
                                'Sig number',
                                #'Wilcoxon Rank Sum',
                                'KS',
                                'Weighted KS'))+
    theme_minimal(base_size = 20)+
    theme(axis.text.x = element_text(vjust=10),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust =title_adj))
} else {
  p =summary2plot %>% 
    filter(!null_type=='random_brain') %>% 
    mutate(null_type_sort=factor(null_type, level = c('random_gene',
                                                      'spin_brain'))) %>% 
    filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
    mutate(stat_sort=factor(stat, level=c('mean',
                                          'meanabs',
                                          'meansqr',
                                          'maxmean',
                                          'median',
                                          #'sign_test',
                                          'sig_n',
                                          #'rank_sum',
                                          'ks_orig',
                                          'ks_weighted'))) %>% 
    mutate(vlineposition=as.numeric(factor(null_type, levels = c('random_gene',
                                                                 'random_brain',
                                                                 'spin_brain'))),
           vlineposition=vlineposition+.5,
           vlineposition=ifelse(vlineposition == max(vlineposition), NA, vlineposition)) %>% 
    ggplot(aes(x=null_type_sort,y=mean_val,fill=stat_sort))+
    geom_bar(stat="identity",position="dodge2") +
    geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit),position="dodge2")+
    ylab(ylab2show)+
    xlab('Null model type')+
    geom_hline(aes(yintercept=0.05),linetype="dashed")+
    ggtitle(title2show)+
    geom_vline(aes(xintercept = vlineposition),color='grey',alpha=0.5)+
    scale_x_discrete(labels= c('Competetive null models',
                               'Self-contained null models'))+
    scale_fill_viridis(discrete = TRUE,
                       name='Test statistic',
                       labels=c('Mean',
                                'Meanabs',
                                'Meansqr',
                                'Maxmean',
                                'Median',
                                #'Sign Test',
                                'Sig Number',
                                #'Wilcoxon Rank Sum',
                                'KS',
                                'Weighted KS'))+
    theme_minimal(base_size = 20)+
    theme(axis.text.x = element_text(vjust=10),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust =title_adj))
  }   
  return(p)
}



# get gene set information: co-expression and n
report_geneSetList <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                               type=c('RandomSet500_int', # intersection with abagen_r04
                                       'GO_BP_int', # intersection with abagen_r04
                                       'GO_CC_int', # intersection with abagen_r04
                                       'GO_MF_int', # intersection with abagen_r04
                                       'SynGO_int',# intersection with abagen_r04
                                       'Celltype_Shin')){
type=match.arg(type)
geneSetList=get_geneSetList(data_path,type=type)
df1=data.frame(genes=unlist(lapply(geneSetList, paste0,collapse='/')))
df2=data.frame(size=unlist(lapply(geneSetList,length)))
df=cbind(df2,df1)
return(df)
}

get_geneSetList_info <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                                  geneSet_type=c('RandomSet500_int', # intersection with abagen_r04
                                         'GO_BP_int', # intersection with abagen_r04
                                         'GO_CC_int', # intersection with abagen_r04
                                         'GO_MF_int', # intersection with abagen_r04
                                         'SynGO_int',# intersection with abagen_r04
                                         'Celltype_Shin'),
                                  geneData_type=c('abagen_r02',
                                         'abagen_r04',
                                         'abagen_r06',
                                         'shin2stage'),
                                  geneData_hem=c('L','R','B')){
  
  # geneSet_type='RandomSet500_int'
  # geneData_type='abagen_r04'
  # geneData_hem='L'
 geneSet_type=match.arg(geneSet_type)
 geneData_type=match.arg(geneData_type)
 geneData_hem=match.arg(geneData_hem)
 geneSetList=get_geneSetList(data_path=data_path,type=geneSet_type)
 gene_data=load_dk_gene_data(data_path=data_path, type=geneData_type,hem=geneData_hem)
 # co-expression
 coexp.mean=lapply(geneSetList, function(x){
            coexp.matrix=cor(gene_data[,x],method='pearson')
            coexp.matrix.upper=coexp.matrix[upper.tri(coexp.matrix, diag = FALSE)]
            coexp.mean=mean(coexp.matrix.upper)
            return(coexp.mean)})
df1=data.frame(coexp_mean=unlist(coexp.mean))
df2=data.frame(size=unlist(lapply(geneSetList,length)))
df_coexp=cbind(df2,df1)
return(df_coexp)
}

#info=get_geneSetList_coexp(geneData_type = 'abagen_r04') 


# get Moran's I and dip test for a set of brain
get_brain_info <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                           brainData_type=c('no_spatial',
                                  'random_spatial003', 
                                  'random_spatial002', 
                                  'random_spatial001',
                                  'matched_spatial003',
                                  'real_brain'),
                           geneData_type=c('abagen_r02',
                                           'abagen_r04',
                                           'abagen_r06',
                                           'shin2stage'),
                           geneData_hem=c('L','R','B')){
  # brainData_type='random_spatial003'
  # geneData_type='abagen_r04'
  # geneData_hem='L'
  
  brainData_type=match.arg(brainData_type)
  geneData_type=match.arg(geneData_type)
  geneData_hem=match.arg(geneData_hem)
  gene_data=load_dk_gene_data(data_path=data_path, type=geneData_type,hem=geneData_hem)
  brain_data=get_brain_data(data_path=data_path,type=brainData_type,col_idx = NA) #NA means all cols
  
  # get Moran's I index
  df.moran=get_column_moran_idx(brain_data)
  # get dip test
  geneList=corr_brain_gene(gene_data=gene_data,
                           brain_data=brain_data)
  # modetest_stat=apply(geneList, 2, function(x) {tmp=multimode::modetest(x,method='ACR')
  #                                                 tmp$statistic})
  modetest_stat=apply(geneList, 2, diptest::dip) # dip is quicker and other approach did not give different subsequent results 
  df.dip=data.frame(modetest_stat=modetest_stat)
  
  # get pval and ifsig
  modetest_pval=apply(geneList, 2, get_dip_pval)
  df.dip_pval=data.frame(modetest_pval=modetest_pval)
  df.dip_ifsig=data.frame(modetest_ifsig=as.factor(as.numeric(modetest_pval<0.05)))
  
  # pos neg mean dist
  pos_neg_dist=apply(geneList, 2, get_pos_neg_dist)
  df.pos_neg_dist=data.frame(pos_neg_dist=pos_neg_dist)
  
  df_brain_info=cbind(df.moran,
                      df.dip,
                      df.dip_pval,
                      df.dip_ifsig,
                      df.pos_neg_dist)
  rownames(df_brain_info)=stringr::str_replace(rownames(df_brain_info),'null_','brain_') # remove null_ to match psig_df
  return(df_brain_info)
}



get_dip_pval <- function(x){
  dip_result=diptest::dip.test(x, simulate.p.value = F)
  return(dip_result$p.value)
}

get_pos_neg_dist <- function(x,
                             by=c('mode','mean')){
  by=match.arg(by)
  if (by=='mean'){
  pos_mean=mean(x[x>=0])
  neg_mean=mean(x[x<0])
  pos_neg_dist=pos_mean-neg_mean
  } else if (by=='mode'){
  pos_dense=density(x[x>=0])
  neg_dense=density(x[x<0])
  pos_mo = pos_dense$x[which.max(pos_dense$y)]
  neg_mo = neg_dense$x[which.max(neg_dense$y)]
  pos_neg_dist=pos_mo-neg_mo
  }
return(pos_neg_dist)
}


correlate_psig_with_info <- function(psig_df,
                                     info,
                                     var2test=c('coexp_mean',
                                                'moran_idx',
                                                'modetest_stat',
                                                'modetest_pval',
                                                'modetest_ifsig',
                                                'pos_neg_dist')){
  #====test nest_by=geneSet=========
  # psig_df=psig.list[[2]]
  # info=get_geneSetList_coexp(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
  # geneSet_type='RandomSet500_int',
  # geneData_type='abagen_r04')
  # var2test='coexp_mean'
  
  #====test nest_by=brain=========
  # psig_df=psig.list[[1]]
  # info=get_brain_info(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
  #                     brainData_type = 'random_spatial003',
  #                     geneData_type = 'abagen_r04')
  # var2test='modetest_stat'
  # var2test='moran_idx'
  # var2test='pos_neg_dist'

  var2test=match.arg(var2test)
  if (!var2test %in% colnames(info)){stop('vart2test is not in the info')}
  psig_df=psig_df %>% tibble::column_to_rownames(colnames(psig_df)[1])
  #if (!identical(order(rownames(psig_df)),order(rownames(info)))){stop('psig_df and info are not matched')}
  
  if (var2test=='coexp_mean'){ # keep size infomation for geneSet
  info2add= info %>%  # here another varaible is ignore (e.g., size for coexp_info)
        rename(var2test=var2test) %>% 
        select(var2test,size)
  res.nested=merge(psig_df,info2add,by='row.names') %>% 
            tibble::column_to_rownames('Row.names') %>% 
            pivot_longer(-c(size,var2test),names_to = 'stat',values_to = 'psig') }
  else {
    info2add= info %>%  
      rename(var2test=var2test) %>% 
      select(var2test)
    res.nested=merge(psig_df,info2add,by='row.names') %>% 
      tibble::column_to_rownames('Row.names') %>% 
      pivot_longer(-c(var2test),names_to = 'stat',values_to = 'psig')}
  
  #res.nested=merge(psig_df,info2add,by='row.names') %>% select(meanabs,var2test)
  
   res.nested=res.nested %>% 
            group_by(stat) %>% 
            nest() %>% 
            mutate(lm_model=map(data,~lm(psig~var2test,.x)),
                   lm_tidy=map(lm_model, tidy),
                   lm_glance=map(lm_model, glance))
  return(res.nested)
}



report_res.nested <- function(res.nested){
    res.df=res.nested %>% 
           unnest(lm_glance) %>% 
           select(-c(statistic, p.value)) %>% # remove duplicated variables// used those in lm_tidy
          filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
           mutate(stat_title=recode(stat,mean='Mean',
                               median='Median',
                               meanabs='Meanabs',
                               meansqr='Meansqr',
                               maxmean='Maxmean',
                               #sig_n='Sig Number',
                               sign_test='Sign Test', 
                               #rank_sum='Wilcoxon Rank Sum',
                               ks_orig='KS',
                               ks_weighted='Weighted KS')) %>% 
          unnest(lm_tidy) %>% filter(term=='var2test'|term=='var2test1')  %>% 
          ungroup() %>% 
          select(-c(data,lm_model)) %>% 
          select(stat_title,statistic,p.value,r.squared) %>% 
          mutate(r.squared=scales::percent(r.squared,accuracy=0.01)) %>% 
          rename('Test statistic'=stat_title,
                 't value'=statistic,
                 'p value'=p.value,
                 'R-squared'=r.squared)
    return(res.df)
}


plot_res.nested <- function(res.nested,
                            annot_position=c(-0.072, 0.55),
                            xlim2show=c(-0.1,0.15),
                            ylim2show=c(-0.05,0.55)){
    res.nested=res.nested %>% 
      filter(!(stat %in% c('sign_test','rank_sum'))) %>% 
      mutate(stat_title=recode(stat,mean='Mean',
                                               median='Median',
                                               meanabs='Meanabs',
                                               meansqr='Meansqr',
                                               maxmean='Maxmean',
                                               sig_n='Sig Number',
                                               #sign_test='Sign Test', 
                                               #rank_sum='Wilcoxon Rank Sum',
                                               ks_orig='KS',
                                               ks_weighted='Weighted KS'),
                              rsqr=map(lm_glance, function(x) x$r.squared)) %>% 
                              unnest(rsqr) %>% 
                              mutate(rsqr=scales::percent(rsqr,accuracy=0.01)) %>% 
    mutate(scatter=pmap(list(data,stat_title,rsqr),plot_psig_info,annot_position=annot_position,xlim2show=xlim2show,ylim2show=ylim2show))
   plot.list=res.nested$scatter
return(plot.list)   
}



plot_psig_info <- function(df2plot,
                           title2show,
                           rsqr2show, 
                           annot_position=c(-0.072, 0.55),
                           xlim2show=c(-0.1,0.15),
                           ylim2show=c(-0.05,0.55)){
  if ('size' %in% colnames(df2plot)){
  p=ggplot(df2plot,aes(var2test, psig, colour=as.factor(size)))+
    geom_point(position = position_jitter(w = 0.01, h = 0.01),
               alpha=0.6,
               size=2)+
    scale_color_viridis(discrete = TRUE)+
    geom_hline(yintercept=0.05,linetype="dashed", color = "darkgrey",size=0.8) +
    ylim(ylim2show)+
    xlim(xlim2show)+
    ggtitle(title2show)+
    annotate(geom='text',x=annot_position[1], y=annot_position[2], label=sprintf('R-sqaured: %s',rsqr2show),color="black",size=4)+
    theme_minimal(base_size = 12)+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5))
  }else{
    p=ggplot(df2plot,aes(var2test, psig))+
      geom_point(position = position_jitter(w = 0.01, h = 0.01),
                 alpha=0.6,
                 size=2,
                 color='#21918c')+
      geom_hline(yintercept=0.05,linetype="dashed", color = "darkgrey",size=0.8) +
      ylim(ylim2show)+
      xlim(xlim2show)+
      ggtitle(title2show)+
      annotate(geom='text',x=annot_position[1], y=annot_position[2], label=sprintf('R-sqaured: %s',rsqr2show),color="black",size=4)+
      theme_minimal(base_size = 12)+
      theme(legend.position = "none",
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5))
  }
  
  return(p)  
}


















