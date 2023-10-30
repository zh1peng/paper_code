args = commandArgs(trailingOnly=TRUE)
slrum_idx=as.numeric(args[1]) # 1,2,3,4,5,6,7,8,9,10
atlas=args[2] # desikan, schaefer100, schaefer200
rdonor=args[3] # r0.2, r0.4, r0.6
gs_type=args[4] # MF, Sim, SynGO
brain_type=args[5] # sample_data, sim_nospatial, sim_spatial0.03, sim_spatial0.02, sim_spatial0.01, real_brain
null_type=args[6] # spin_brain, random_gene, random_gene_coexp, random_gene_subset, spin_random_mixed
cor_type=args[7] # pearson, spearman, loo, pls1, pls1w


# print all input parameters
cat('==========\n')
cat('input parameters:\n')
cat(sprintf('atlas: %s\n',atlas))
cat(sprintf('rdonor: %s\n',rdonor))
cat(sprintf('gs_type: %s\n',gs_type))
cat(sprintf('brain_type: %s\n',brain_type))
cat(sprintf('null_type: %s\n',null_type))
cat(sprintf('cor_type: %s\n',cor_type))
cat(sprintf('slrum_idx: %s\n',slrum_idx))
cat('==========\n')

# setup paths 
data_path='/gpfs1/home/z/c/zcao4/revision_code/data'
code_path='/gpfs1/home/z/c/zcao4/revision_code/functions'
res_path='/gpfs1/home/z/c/zcao4/revision_code/results'
# source functions
source(sprintf( '%s/data_functions.R',code_path))
source(sprintf( '%s/spinning_brain_functions.R',code_path))
source(sprintf( '%s/resampling_gene_functions.R',code_path))
source(sprintf( '%s/cor_functions.R',code_path))
source(sprintf( '%s/vacc_function.R',code_path))

# create path
save_path=file.path(res_path,sprintf('%s_%s_%s_%s_%s_%s',atlas,rdonor,brain_type,gs_type,null_type,cor_type))
dir.create(save_path, showWarnings = FALSE)

savename=sprintf('%s/sim_%s_res.csv',save_path, slrum_idx)
if (file.exists(savename)){
  stop('file exists')
}



tt=Sys.time()
if (null_type %in% c('random_gene_coexp', 'random_gene_subset')){

  constrain=switch(null_type,
                   'random_gene_coexp'='match_coexp',
                   'random_gene_subset'='gene_subset')
sampled_geneSetList=load_sampled_GeneSets(sprintf('%s/sampled_GeneSets',data_path),
                      atlas=atlas,
                      rdonor=rdonor,
                      gs_type=gs_type,
                      constrain=constrain)
} else {
    sampled_geneSetList=NULL
}



df=vacc_function(data_path=data_path,
                 atlas=atlas,
                 rdonor=rdonor,
                 gs_type=gs_type,
                 brain_type=brain_type,
                 null_type=null_type,
                 cor_type=cor_type,
                 sampled_geneSetList=sampled_geneSetList,
                 slrum_idx=slrum_idx)
write.csv(df,savename,row.names = F)