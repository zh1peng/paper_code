library(ape)
# simulate brain data with desired Moran range
generate_data_with_moran <- function(atlas = c('desikan',
                                              'schaefer100',
                                              'schaefer200',
                                              'schaefer300'),
                                     moran.m = 0.03,
                                     moran.sd = 0.001,
                                     n.pool = 100000,
                                     n.sample = 1000,
                                     method = c('runif', 'rnorm'),
                                     seed2use = 2022) {
  atlas=match.arg(atlas)
  method=match.arg(method)
  set.seed(seed2use)
  
  n.region <- switch(atlas,
    'desikan'       = 34,
    'schaefer100'   = 50,
    'schaefer200'   = 100,
    'schaefer300'   = 150,
    NA_integer_
  )

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
  
  random_brain_moran = get_column_moran_idx(weight_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights', 
                                            atlas=atlas,
                                            df=random_brain_data)
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
  colnames(sampled_random_brain_data) = paste0('sim_', c(1:n.sample))
  return(sampled_random_brain_data)
}


get_column_moran_idx <- function(weight_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/perm_id_weights', 
                                 atlas=c('desikan',
                                          'schaefer100',
                                          'schaefer200',
                                          'schaefer300'),
                                 df
                                  ){
  atlas=match.arg(atlas)
  weight_file=sprintf('%s/%s_Moran_w.rds',weight_path,atlas)
  #check if weight_file exist
    if (!file.exists(weight_file)) {
        stop(sprintf('weight_file %s does not exist',weight_file))
    }
  w=readRDS(weight_file)
  moran_idx=apply(df, 2, function(x){
                            temp=Moran.I(x,w)
                            return(temp$observed)
                            }) 
  return(as.data.frame(moran_idx))
}



generate_random_data <- function(atlas = c('desikan',
                                            'schaefer100',
                                            'schaefer200',
                                            'schaefer300'),
                                            n.sample = 1000,
                                            seed2use = 2022){


  atlas=match.arg(atlas)
  set.seed(seed2use)
  
  n.region <- switch(atlas,
    'desikan'       = 34,
    'schaefer100'   = 50,
    'schaefer200'   = 100,
    'schaefer300'   = 150,
    NA_integer_
  )                                              
random_brain = matrix(runif(n.region * n.sample, min = -1, max = 1), nrow =n.region)
colnames(random_brain)=paste0('null_', c(1:n.sample))
return(random_brain)
}

get_region_labels<-function(label_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainInfo/labels/',
                            atlas=c('desikan',
                                    'schaefer100',
                                    'schaefer200',
                                    'schaefer300')){
    atlas=match.arg(atlas)
    label_file=sprintf('%s/%s_labels.csv',label_path,atlas)
    # check if label_file exist
    if (!file.exists(label_file)) {
        stop(sprintf('label_file %s does not exist',label_file))
    }

    hem='L'
     if (atlas=='desikan'){
    search_pattern=sprintf('^%s_',hem)
    }else if (atlas %in% c('schaefer100','schaefer200','schaefer300')){
        search_pattern=sprintf('%sH', hem)
    }

    label2add=read.csv(label_file) %>% filter(grepl(search_pattern,label))
    return(label2add)
}

# # generate brain data
# for (atlas_i in c('desikan')){
#     for (moran_idx in c(NA,0.03, 0.02,0.01)){
#         if (is.na(moran_idx)){
#             random_brain=generate_random_data(atlas=atlas_i)
#             savename=sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainDat/%s_sim_nospatial.csv',atlas_i)
#         } else {
#             random_brain=generate_data_with_moran(atlas=atlas_i,moran.m=moran_idx)
#             savename=sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainDat/%s_sim_spatial%s.csv',atlas_i,moran_idx)
#         }
#         label2add=get_region_labels(atlas=atlas_i)
#         rownames(random_brain)=label2add$label
#         write.csv(random_brain,savename)
#     }
# }
# atlas_i='schaefer200'
# colMeans(get_column_moran_idx(atlas='schaefer200',df=random_brain))


# # generate brain data
# for (atlas_i in c('schaefer200')){
#     for (moran_idx in c(0.03)){
#         if (is.na(moran_idx)){
#             random_brain=generate_random_data(atlas=atlas_i)
#             savename=sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainDat/%s_sim_nospatial.csv',atlas_i)
#         } else {
#             random_brain=generate_data_with_moran(atlas=atlas_i,moran.m=moran_idx,n.pool = 500000,seed2use = 2023)
#             savename=sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/BrainDat/%s_sim_spatial%s.csv',atlas_i,moran_idx)
#         }
#         label2add=get_region_labels(atlas=atlas_i)
#         rownames(random_brain)=label2add$label
#         write.csv(random_brain,savename)
#     }
# }
#atlas_i='schaefer200'
#colMeans(get_column_moran_idx(atlas='schaefer200',df=random_brain))

