library(ape)
library(spdep)

# assuming code_path is set 
# code_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code'
w=readRDS(file.path(code_path,'/data/Moran_w.rds')

get_column_moran_idx <- function(df){
  moran_idx=apply(df, 2, get_moran_idx_ape, weight=.GlobalEnv$w) # read w first
  return(as.data.frame(moran_idx))
}
get_moran_idx_ape <- function(x,weight){
  temp=Moran.I(x,weight)
  return(temp$observed)
}

normalize_weight <- function(weight){
  ROWSUM = rowSums(weight)
  ROWSUM[ROWSUM == 0] = 1
  weight.norm = weight/ROWSUM
  return(weight.norm)
}


get_column_localmoran <- function(df,stat=c('Ii','Pr(z != E(Ii))')){
  stat=match.arg(stat)
  moran.local=apply(df,2,function(x) {temp=get_localmoran_spdep(x,weight=.GlobalEnv$w)
                                     return(temp[,stat])})
  return(as.data.frame(moran.local))
}

get_localmoran_spdep <- function(x,weight){
  weight=normalize_weight(weight)
  listw=mat2listw(weight)
  temp=localmoran(x,listw)
  return(temp)
}


get_column_gearyC <- function(df){
  gearyC=apply(df,2,function(x) {temp=get_gearyC_spdep(x,weight=.GlobalEnv$w)
                                        return(temp$estimate[1])})
  return(as.data.frame(gearyC))
}

get_gearyC_spdep <- function(x,weight){
  weight=normalize_weight(weight)
  listw=mat2listw(weight)
  temp=geary.test(x,listw)
  return(temp)
}

# find path of nearest region
#data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data'
#lh_centroid_data=read.csv(sprintf('%s/lh_centroid.csv',data_path))
#centroid.l=lh_centroid_data %>%select(starts_with('centroid'))
#rownames(centroid.l)=lh_centroid_data$Row
#orig.label=read.csv(sprintf('%s/dk_standard_var_order.csv',data_path)) %>% filter(grepl('^L_',label))
#orig.label=orig.label$label

get_reorder_idx <- function(centroid,
                            start='L_frontalpole',
                            orig.label,
                            method=c('nearest_path','nearest_region','overall_dist')){
  method=match.arg(method)
  
  if (method == 'nearest_path') {
    nearest_region = start
    for (iter in c(1:dim(centroid)[1] - 1)) {
      dist.mat = as.matrix(dist(centroid))
      rowdata = sort(dist.mat[nearest_region[iter],])
      found_region = names(rowdata)[!(names(rowdata) %in% nearest_region)][1]
      nearest_region = c(nearest_region, found_region)
    }
    reorder.idx = match(orig.label, nearest_region)
  } else if (method == 'nearest_region'){
    dist.mat = as.matrix(dist(centroid))
    rowdata = sort(dist.mat[start, ])
    reorder.idx = match(orig.label,names(rowdata))
  } else if (method == 'overall_dist'){
    dist.mat = as.matrix(dist(centroid))
    rowdata=sort(rowSums(dist.mat))
    reorder.idx = match(orig.label,names(rowdata))
  }
  return(reorder.idx)
}

# reorder.idx=get_reorder_idx(centroid=centroid.l,
#                 start ='L_frontalpole',
#                 orig.label,
#                 method='nearest_path')


            
# alpha=0.3
# df=read.csv(sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/new_code/random_brain_maps/alpha%s_lh_simbrain_DK.csv',alpha),
#             row.names = 1)
# 
# df.moran=get_df_moran_idx(df)    
# hist(df.moran$moran_idx,breaks = 100)

# approach 1 select brain maps with desired range
generate_data_with_moran <- function(n.region = 34,
                                     moran.m = 0.03,
                                     moran.sd = 0.001,
                                     n.pool = 500000,
                                     n.sample = 1000,
                                     method = c('runif', 'rnorm'),
                                     seed2use = 2022) {
  method=match.arg(method)
  set.seed(seed2use)
  
  
  if (moran.m < 0) {
    stop('moran.m<0 is not supported')
  }
  
  if (method == 'runif') {
    random_brain_data = matrix(runif(n.region * n.pool, min = -1, max = 1), nrow =n.region)
  } else if (method == 'rnorm') {
    random_brain_data = matrix(rnorm(n.region* n.pool, mean=0, sd = 1), nrow =n.region)
  }
  
  random_brain_data = as.data.frame(random_brain_data)
  colnames(random_brain_data) = paste0('null_', c(1:n.pool))
  
  random_brain_moran = get_column_moran_idx(random_brain_data)
  random_brain_moran.filtered = random_brain_moran[random_brain_moran$moran_idx >0, , drop = F]
  random_brain_moran.filtered$prob = dnorm(random_brain_moran.filtered$moran_idx,
                                           mean = moran.m,
                                           sd = moran.sd)
  sampled_moran = sample(random_brain_moran.filtered$moran_idx,
                         n.sample,
                         prob = random_brain_moran.filtered$prob)
  
  # find their data
  
  sampled_cols = rownames(random_brain_moran.filtered)[random_brain_moran.filtered$moran_idx %in% sampled_moran]
  sampled_random_brain_data = random_brain_data[, sampled_cols]
  
  # check moran idx
  #hist(get_column_moran_idx(sampled_random_brain_data)$moran_idx)
  
  #rename colnames
  colnames(sampled_random_brain_data) = paste0('null_', c(1:n.sample))
  
  return(sampled_random_brain_data)
}

## approach 2 find shuffled data to match the desired Moran value
generate_data_with_moran1 <- function(moran.m=0.03,
                                      moran.sd=0.001,
                                      n.region=34,
                                      n.sample=1000,
                                      stop.iter=5000,
                                      weight,
                                      orig=NULL){
  
  moran.upper = moran.m + moran.sd
  moran.lower = moran.m - moran.sd
  res = matrix(NA, ncol = n.sample, nrow = n.region)
  orig= matrix(NA, ncol = n.sample, nrow = n.region)
  found.n = 0
  while (found.n < n.sample ) {
    mI = 0
    iter=0
    x = runif(34,-1, 1)
    while (mI > moran.upper |
           mI < moran.lower | 
           iter < stop.iter)
    {
      iter=iter+1
      x.shuffle = x[sample(34)]
      mI = get_moran_idx_ape(x.shuffle, weight)
    }
    
    found.n = found.n + 1
    orig[,found.n]=x
    res[, found.n] = x.shuffle
    cat(sprintf('found: %s mI: %s with %s iteration \n',found.n,mI,iter))
  }
  return(list(orig.x=orig,
              res.x=res))
}

# approach 2 but with a simple function
match_moran <- function(x,
                        weight,
                        moran.m=0.03,
                        moran.sd=0.001){
  # example:
  # x=runif(34,-1,1)
  # match_moran(x,weight=w)
  
  mI=0
  iter=0
  moran.upper = moran.m + moran.sd
  moran.lower = moran.m - moran.sd
  
  while (mI > moran.upper |
         mI < moran.lower | 
         iter < stop.iter)
  {
    iter=iter+1
    x.shuffle = x[sample(34)]
    mI = get_moran_idx_ape(x.shuffle, weight)
  }
  cat(sprintf('found mI: %s with %s iterations \n',mI,iter))
  return(x.shuffle)
}

            
            