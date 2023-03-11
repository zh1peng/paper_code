# vacc function

args = commandArgs(trailingOnly=TRUE)
slrum.array.idx=as.numeric(args[1])
random_brain_type=args[2] # random_spatial003, matched_spatial003, no_spatial, random_spatial002, random_spatial001
gene_data_type= args[3] # abagen_r02,abagen_r04, abagen_r06, shin2stage
gene_set_type= args[4] #  'RandomSet_test  'RandomSet500', 'GO_BP', 'SynGO'
null_model_type=args[5] # 'spin_brain' 'random_brain' 'random_gene'
null_model_n=as.numeric(args[6]) # 5000 or 10000

data.path='/gpfs1/home/z/c/zcao4/geneset_sim_analysis/data'
res.path='/gpfs1/home/z/c/zcao4/geneset_sim_analysis/all_sim'



if(! random_brain_type %in% c('random_spatial003', 
                              'matched_spatial003', 
                              'no_spatial', 
                              'random_spatial002', 
                              'random_spatial001',
                              'real_brain')){
  stop('unknown brain data')
}

if (!gene_set_type %in% c('RandomSet_test', # generate 10 random set based gene_data
                          'RandomSet500', # generate 500 random sets based gene_data
                          'RandomSet500_int', # intersection with abagen_r04
                          'GO_BP_int', # intersection with abagen_r04
                          'GO_CC_int', # intersection with abagen_r04
                          'GO_MF_int', # intersection with abagen_r04
                          'SynGO_int',# intersection with abagen_r04
                          'Celltype_Shin')) {
  stop('unknown geneSetList')
}

if (!gene_data_type %in% c('abagen_r02', 
                           'abagen_r04', 
                           'abagen_r06', 
                           'shin2stage')) {
  stop('unknown gene type')
}

if (!null_model_type %in% c('spin_brain', 
                            'random_brain', 
                            'random_gene')) {
  stop('unknown null model type')
}

cat(sprintf('==========\n%s\n%s\n%s\n%s\n%s\n',gene_set_type, gene_data_type,
            random_brain_type, null_model_type, null_model_n))
save.path=file.path(res.path,sprintf('%s_%s_%s_%s_%s',gene_set_type, gene_data_type,
                                     random_brain_type, null_model_type, null_model_n))

if (file.exists(sprintf('%s/sim_%s_res.csv',save.path, slrum.array.idx))){
  stop('file exists')
}

dir.create(save.path, showWarnings = FALSE)




source(file.path(data.path,'cor_functions.R'))
source(file.path(data.path,'data_functions.R'))

brain_data=get_brain_data(data_path=data.path,
                         type=random_brain_type,
                         col_idx=slrum.array.idx)

gene_data=load_dk_gene_data(data_path=data.path,
                            type=gene_data_type)

geneSetList=get_geneSetList(data.path,
                            gene_set_type)
# if (gene_set_type=='RandomSet_test'){
#   geneSetList=simulate_geneSetList(gene_data,sim.rep = 1)
# } else if (gene_set_type=='RandomSet500'){
#   geneSetList=simulate_geneSetList(gene_data)
# } else if (gene_set_type=='RandomSet500_int'){
#   geneSetList=readRDS(sprintf('%s/RandomSet500_int.rds',data.path))
# } else if (gene_set_type=='GO_BP_int'){
#   geneSetList=readRDS(sprintf('%s/GO_BP_20_200_2247_int.rds',data.path))
# }else if (gene_set_type=='GO_CC_int'){
#   geneSetList=readRDS(sprintf('%s/GO_CC_20_200_305_int.rds',data.path))
# }else if (gene_set_type=='GO_MF_int'){
#   geneSetList=readRDS(sprintf('%s/GO_MF_20_200_316_int.rds',data.path))
# }else if (gene_set_type=='SynGO_int'){
#   geneSetList=readRDS(sprintf('%s/SynGO_20_200_49_int.rds',data.path))
# }else if (gene_set_type=='Celltype_Shin'){
#   geneSetList=readRDS(sprintf('%s/Celltype_shin_9.rds',data.path))
# }



geneList.true=corr_brain_gene(gene_data,brain_data)
if (null_model_type=='random_gene') {
  geneList.null = resampling_geneList(geneList.true,
                                      perm.n = null_model_n)
} else {
  perm.id = load_34dk_permid(data_path = data.path,
                             type = null_model_type,
                             perm.n = null_model_n)
  null_brain_data = generate_null_brain_data(brain_data, perm.id)
  geneList.null = corr_brain_gene(gene_data, null_brain_data)
}

cat('==========\n')
cat('start testing \n')
tt=Sys.time()
res=list()
for (method2test in c('mean','median','meanabs','meansqr','maxmean','sig_n',
                      'sign_test','rank_sum','ks_orig','ks_weighted')){
  cat(paste0(method2test,'....\n'))
  score.true=aggregate_geneSetList(geneSetList, geneList.true,method=method2test)
  score.null=aggregate_geneSetList(geneSetList, geneList.null,method=method2test)
  score.null.upper25=lapply(score.null, function(x){as.numeric(quantile(x,prob=0.975))})
  score.null.lower25=lapply(score.null, function(x){as.numeric(quantile(x,prob=0.025))})
  
  score.null.upper5=lapply(score.null, function(x){as.numeric(quantile(x,prob=0.95))})
  score.null.lower5=lapply(score.null, function(x){as.numeric(quantile(x,prob=0.05))})
  
  pvals=caculate_pvals(score.true,score.null)
  res[[method2test]][['true']]=score.true
  res[[method2test]][['null_upper25']]=score.null.upper25
  res[[method2test]][['null_lower25']]=score.null.lower25
  res[[method2test]][['null_upper5']]=score.null.upper5
  res[[method2test]][['null_lower5']]=score.null.lower5
  res[[method2test]][['pvals']]=pvals
  rm(score.true,score.null,pvals,score.null.upper25,score.null.lower25,score.null.upper5,score.null.lower5)
}


df.list=list()
for (var2extract in c('true','null_upper25','null_lower25','null_upper5','null_lower5','pvals')){
  df.list[[var2extract]]= as.data.frame(lapply(res, function(x){do.call(rbind,x[[var2extract]])}))
}
df2save=do.call(cbind,df.list) %>%  tibble::rownames_to_column("geneSet") %>% mutate(null_brain=slrum.array.idx)
write.csv(df2save,sprintf('%s/sim_%s_res.csv',save.path, slrum.array.idx),row.names = F)
Sys.time()-tt

# saveRDS(res,file = sprintf('%s/sim_%s.rds',save.path, slrum.array.idx))



