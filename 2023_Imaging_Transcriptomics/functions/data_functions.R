# load data functions
library(dplyr)
load_dk_gene_data <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                              type=c('abagen_r02',
                                     'abagen_r04',
                                     'abagen_r06',
                                     'shin2stage'),
                              hem=c('L','R','B')){
  type=match.arg(type)
  hem=match.arg(hem)
  if ((!type=='shin2stage')&(!hem=='L')){warning('Two right hemspheric regions with missing values are removed')}
  
  csvname=switch (type,
    'abagen_r02'='allgenes_stable_r0.2.csv',
    'abagen_r04'='allgenes_stable_r0.4.csv',
    'abagen_r06'='allgenes_stable_r0.6.csv',
    'shin2stage'='shin2stage.csv')
    gene.df=read.csv(file.path(data_path,csvname), stringsAsFactors = F)

  if (hem=='L'|hem=='R'){
    gene.df=gene.df %>% filter(grepl(sprintf('^%s_',hem),Region))
  }                
  gene.df=gene.df %>% filter(complete.cases(.)) %>% tibble::column_to_rownames('Region')
  gene.mx=as.matrix(gene.df)
  return(gene.mx)
}

# read brain data
load_dk_brain_data <- function(file_name='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data/sample_brain_phenotype.csv',
                               hem=c('L','R','B')){
  if (missing(hem)){hem='L'}
  
  brain.df=read.csv(file_name,stringsAsFactors = F)
  if (hem=='L'|hem=='R'){
    brain.df=brain.df%>% filter(grepl(sprintf('^%s_',hem),Region))
  }             
  brain.df=brain.df %>% tibble::column_to_rownames(var = "Region") 
  colnames(brain.df)='phenotype'
  brain.mx=as.matrix(brain.df)
  return (brain.mx)
}

get_brain_data <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                                  type=c('no_spatial',
                                         'random_spatial003', 
                                          'random_spatial002', 
                                          'random_spatial001',
                                         'matched_spatial003',
                                         'real_brain'),
                                  col_idx=1){
  type=match.arg(type)
  csvname=switch (type,
          'no_spatial' = 'random_brain_nospatial.csv',
          'random_spatial003'='random_brain_spatial003.csv',
          'random_spatial002'='random_brain_spatial002.csv',
          'random_spatial001'='random_brain_spatial001.csv',
          'real_brain'='real_brain.csv')
  df=read.csv(file.path(data_path,csvname),row.names = 1)
  df.mx=as.matrix(df)
  if (is.na(col_idx[1])){col_idx=c(1:dim(df.mx)[2])}
  df.mx=df.mx[,col_idx,drop=F]
  if (length(col_idx)==1){
  colnames(df.mx)='phenotype'
  }
  return(df.mx)
}


get_geneSetList <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                            type=c('RandomSet_test', # generate 10 random set based gene_data
                                   'RandomSet500', # generate 500 random sets based gene_data
                                   'RandomSet500_int', # intersection with abagen_r04
                                   'GO_BP_int', # intersection with abagen_r04
                                   'GO_CC_int', # intersection with abagen_r04
                                   'GO_MF_int', # intersection with abagen_r04
                                   'SynGO_int',# intersection with abagen_r04
                                   'Celltype_Shin')){
type=match.arg(type)
  if (type=='RandomSet_test'){
    geneSetList=simulate_geneSetList(gene_data,sim.rep = 1)
  } else if (type=='RandomSet500'){
    geneSetList=simulate_geneSetList(gene_data)
  } else if (type=='RandomSet500_int'){
    geneSetList=readRDS(sprintf('%s/RandomSet500_int.rds',data_path))
  } else if (type=='GO_BP_int'){
    geneSetList=readRDS(sprintf('%s/GO_BP_20_200_2247_int.rds',data_path))
  }else if (type=='GO_CC_int'){
    geneSetList=readRDS(sprintf('%s/GO_CC_20_200_305_int.rds',data_path))
  }else if (type=='GO_MF_int'){
    geneSetList=readRDS(sprintf('%s/GO_MF_20_200_316_int.rds',data_path))
  }else if (type=='SynGO_int'){
    geneSetList=readRDS(sprintf('%s/SynGO_20_200_49_int.rds',data_path))
  }else if (type=='Celltype_Shin'){
    geneSetList=readRDS(sprintf('%s/Celltype_shin_9.rds',data_path))
  }
return(geneSetList)
}



load_34dk_permid <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/data',
                             type=c('spin_brain','random_brain'),
                             perm.n=5000){
type=match.arg(type)
if (perm.n>10000){stop('10k perm.id is supported')} 
if (type=='spin_brain'){
 perm.id=readRDS(paste(data_path,'34perm_id_spin10k.rds',sep='/'))   
} else if (type=='random_brain'){
 perm.id=readRDS(paste(data_path,'34perm_id_nospin10k.rds',sep='/')) 
}
  return(perm.id[,1:perm.n])
}


generate_null_brain_data <- function(brain_data,perm.id){
  # before put perm.id in the function 
  # need to make sure no duplicate in perm.id.
  region.n=dim(perm.id)[1]
  perm.n=dim(perm.id)[2]
  if (!dim(brain_data)[1]==region.n){stop('The number of regions in brain data nd perm.id is not matched ')}
  # null_brain_data=matrix(NA,nrow=region.n, ncol=perm.n)
  # for (idx in c(1:perm.n)){
  #   tmp_brain_data=brain_data
  #   null_brain_data[,idx]=brain_data[perm.id[,idx]]
  # }
  null_brain_data=sapply(c(1:perm.n),function(idx){
    brain_data[perm.id[,idx]]
  })
  rownames(null_brain_data)=rownames(brain_data)
  colnames(null_brain_data)=paste0('null_',c(1:perm.n))
  return(null_brain_data)
}



resampling_geneList<- function(geneList.true,
                               perm.n=5000){
  if (perm.n>10000){stop('10k perm.id is supported')} 
  geneList.null= replicate(n = 11000, sample(geneList.true,size=dim(geneList.true)[1],replace=F), simplify = T)
  geneList.null=geneList.null[,!duplicated(t(geneList.null))] # remove dupliated column
  geneList.null=geneList.null[,c(1:perm.n)]
  rownames(geneList.null)=rownames(geneList.true)
  colnames(geneList.null)=paste0('null_',c(1:perm.n))
  attr(geneList.null,'is_fisherz')=attr(geneList.true,'is_fisherz')
  attr(geneList.null,'n.region')=attr(geneList.true,'n.region')
  return(geneList.null)
}





simulate_geneSetList <- function(gene_data,
                                 sim.size = seq(20, 200, by=20),
                                 sim.rep = 50) {
  
  # create simulated geneSetList for each [sim.size]
  # and this is repeated for [sim.rep] times with different seeds
  lapply(sim.size,function(size_i){ 
  set.seed(size_i)
  seeds = sample(100000:999999, size = sim.rep, replace = F)
  names(seeds)=paste0('S',seeds)
  lapply(seeds, function(rep_i) { # alternative: sapply: simplify=F, use.names=T
    set.seed(rep_i)
    sample(colnames(gene_data), size_i, replace = F)
  })}) %>% purrr::flatten()}
  






