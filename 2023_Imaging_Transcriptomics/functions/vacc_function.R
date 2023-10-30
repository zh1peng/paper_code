vacc_function<-function(
                data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data',
                atlas=c('desikan',
                        'schaefer100',
                        'schaefer200'),
                rdonor=c('r0.2',
                        'r0.4',
                        'r0.6'),
                gs_type=c('MF',
                        'Sim',
                        'SynGO'),
                brain_type=c('sample_data',
                            'sim_nospatial',
                            'sim_spatial0.03',
                            'sim_spatial0.02',
                            'sim_spatial0.01',
                            'real_brain'),
                null_type=c('spin_brain','random_gene','random_gene_coexp', 'random_gene_subset','spin_random_mixed'),
                cor_type=c('pearson','spearman','loo','pls1','pls1w'),
                sampled_geneSetList=NULL, # any constrained sampled_geneSetList for random_gene?
                slrum_idx=1){
# print all input parameters
cat('==========\n')
cat('input parameters for vacc fun:\n')
cat(sprintf('data_path: %s\n',data_path))
cat(sprintf('atlas: %s\n',atlas))
cat(sprintf('rdonor: %s\n',rdonor))
cat(sprintf('gs_type: %s\n',gs_type))
cat(sprintf('brain_type: %s\n',brain_type))
cat(sprintf('null_type: %s\n',null_type))
cat(sprintf('cor_type: %s\n',cor_type))
cat(sprintf('is null sampled_geneSetList: %s\n',is.null(sampled_geneSetList)))
cat(sprintf('slrum_idx: %s\n',slrum_idx))
cat('==========\n')

brain_data=load_BrainDat(data_path=sprintf('%s/BrainDat',data_path),
                        atlas=atlas,
                        type=brain_type,
                        col_idx=slrum_idx)
gene_data=load_GeneExp(data_path=sprintf('%s/GeneExp',data_path),
                        atlas=atlas,
                        rdonor=rdonor)
geneSetList=load_GeneSets(data_path=sprintf('%s/GeneSets',data_path),
                        atlas=atlas,
                        rdonor=rdonor,
                        gs_type=gs_type)

geneList.true=corr_brain_gene(gene_data,brain_data,method=cor_type)


if (null_type=='random_gene') {
  geneList.null = resampling_geneList(geneList.true)
} else if (null_type=='spin_brain'){
  perm.id = load_permid(data_path=sprintf('%s/BrainInfo/perm_id_weights',data_path),
                        atlas=atlas,
                        type='spin_brain')
  null_brain_data = generate_null_brain_data(brain_data, perm.id)
  geneList.null = corr_brain_gene(gene_data, null_brain_data,method=cor_type)
} else if (null_type=='spin_random_mixed'){
  perm.id = load_permid(data_path=sprintf('%s/BrainInfo/perm_id_weights',data_path),
                        atlas=atlas,
                        type='spin_brain')
  null_brain_data = generate_null_brain_data(brain_data, perm.id)
  geneList.tmp = corr_brain_gene(gene_data, null_brain_data,method=cor_type)
  geneList.null = apply(geneList.tmp,2,sample)
  row.names(geneList.null)=row.names(geneList.tmp)
  colnames(geneList.null)=paste0('null_',c(1:5000))
  attr(geneList.null,'is_fisherz')=attr(geneList.tmp,'is_fisherz')
  attr(geneList.null,'n.region')=attr(geneList.tmp,'n.region')
}

cat('==========\n')
cat('start testing \n')
res=list()
for (method2test in c('mean','median','meanabs','meansqr','maxmean','sig_n','ks_orig','ks_weighted')){
if (cor_type %in% c('pls1','loo','pls1w')& method2test=='sig_n'){ # skip sig_n for pls1 and loo
    next
}
cat(paste0(method2test,'....\n'))
  score.true=aggregate_geneSetList(geneSetList, geneList.true,method=method2test)
  
  # caculate null stats
  if (null_type %in% c('random_gene','spin_brain','spin_random_mixed')){ # resampling without constrains
  score.null=aggregate_geneSetList(geneSetList, geneList.null,method=method2test)
  } else { # resampling with constrains
  
  if (is.null(sampled_geneSetList)){
    stop('sampled_geneSetList is required for random_gene_coexp and random_gene_subset')
  }
  score.null=aggregate_geneSetList_with_constrain(geneSetList = geneSetList,
                                                   sampled_geneSetList = sampled_geneSetList,
                                                    geneList = geneList.true,
                                                    method=method2test)  
  }
  # caculate p value based on true and null stats
  if (method2test %in% c('ks_orig','ks_weighted')){
    pvals=caculate_pvals(score.true,score.null,method='split_pos_neg')
  } else {
    pvals=caculate_pvals(score.true,score.null,method='standard')
  }
  res[[method2test]][['pvals']]=pvals
  rm(score.true,score.null,pvals)
}


# format res to a data.frame
df.list=list()
for (var2extract in c('pvals')){
  df.list[[var2extract]]= as.data.frame(lapply(res, function(x){do.call(rbind,x[[var2extract]])}))
}
df2save=do.call(cbind,df.list) %>%  tibble::rownames_to_column("geneSet") %>% mutate(sim_brain=slrum_idx)
return(df2save)
}