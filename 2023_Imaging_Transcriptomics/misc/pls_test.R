library(pls)

args = commandArgs(trailingOnly = TRUE)
slrum_idx=as.numeric(args[1]) # 1,2,3,4,5,6,7,8,9,10
atlas=args[2] # desikan, schaefer100, schaefer200
rdonor=args[3] # r0.2, r0.4, r0.6
gs_type=args[4] # MF, Sim, SynGO
brain_type=args[5] # sample_data, sim_nospatial, sim_spatial0.03, sim_spatial0.02, sim_spatial0.01, real_brain
null_type=args[6] # spin_brain, random_gene, random_gene_coexp, random_gene_subset


# slrum_idx=1
# atlas='desikan'
# rdonor='r0.6'
# gs_type='MF'
# brain_type='sim_spatial0.03'
# null_type='spin_brain'


# print all input parameters
cat('==========\n')
cat('input parameters:\n')
cat(sprintf('atlas: %s\n',atlas))
cat(sprintf('rdonor: %s\n',rdonor))
cat(sprintf('gs_type: %s\n',gs_type))
cat(sprintf('brain_type: %s\n',brain_type))
cat(sprintf('null_type: %s\n',null_type))
cat(sprintf('slrum_idx: %s\n',slrum_idx))
cat('==========\n')

# setup paths 
revision_code_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code'             
data_path=sprintf('%s/data',revision_code_path)
code_path=sprintf('%s/functions',revision_code_path)
res_path=sprintf('%s/results',revision_code_path)

source(sprintf( '%s/data_functions.R',code_path))
source(sprintf( '%s/spinning_brain_functions.R',code_path))

get_true_rmse <- function(gene_data,brain_data,gs){
pls.model=plsr(brain_data~gene_data[,gs], validation='CV')
rmse.vals=drop(unlist(RMSEP(pls.model)$val))
rmse.vals=rmse.vals[,-1] # remove intercept column
rmse.true=min(rmse.vals[2,])
ncomp2use=unname(which.min(rmse.vals[2,]))
return(list(rmse.true=rmse.true, ncomp2use=ncomp2use))
}

get_null_rmse_sampled_gs <- function(gene_data, brain_data, sampled_gs, ncomp2use){

null_rmse=sapply(sampled_gs, function(gs_i) {

  pls.model=plsr(brain_data~gene_data[,gs_i], ncomp=ncomp2use, validation='CV')
  rmse.vals=drop(unlist(RMSEP(pls.model)$val))
  rmse.vals=rmse.vals[,-1,drop=F] # remove intercept column
  rmse.null=rmse.vals[2,ncomp2use]
  return(rmse.null)})
return(null_rmse)
}

get_null_rmse_spin_brain <-function(gene_data, spin_brain_data, gs, ncomp2use){

rmse.null=apply(spin_brain_data, 2, function(brain_data_i) {
  pls.model=plsr(brain_data_i~gene_data[,gs], ncomp=ncomp2use, validation='CV')
  rmse.vals=drop(unlist(RMSEP(pls.model)$val))
  rmse.vals=rmse.vals[,-1,drop=F] # remove intercept column
  rmse.null=rmse.vals[2,ncomp2use]
  return(rmse.null)})
}


# create path
save_path=file.path(res_path,sprintf('%s_%s_%s_%s_%s_PLSR',atlas,rdonor,brain_type,gs_type,null_type))
dir.create(save_path, showWarnings = FALSE)

savename=sprintf('%s/sim_%s_res.csv',save_path, slrum_idx)
if (file.exists(savename)){
  stop('file exists')
}


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


if (!identical(rownames(brain_data), rownames(gene_data))) {
    stop('Regions in brain_data and gene_data are not matched')
  }

# prepare neccesary variables for null models
if (null_type %in% c('random_gene_coexp', 'random_gene_subset')){

constrain=switch(null_type,
                   'random_gene_coexp'='match_coexp',
                   'random_gene_subset'='gene_subset')
sampled_geneSetList=load_sampled_GeneSets(sprintf('%s/sampled_GeneSets',data_path),
                      atlas=atlas,
                      rdonor=rdonor,
                      gs_type=gs_type,
                      constrain=constrain)
  if(!identical(names(sampled_geneSetList), names(geneSetList))){
    stop('names of sampled_geneSetList and geneSetList are not matched')
 }

} else if (null_type=='random_gene'){
  set.seed(slrum_idx)
  perm.n=5000
  sampled_geneSetList=lapply(geneSetList,function(gs){
    sampled_gs=lapply(1:perm.n+100, function(x) {sample(colnames(gene_data),size=length(gs),replace=F)})
    sampled_gs=sampled_gs[!duplicated(sampled_gs)]
    sampled_gs=sampled_gs[1:perm.n]# remove duplciated sampled_gs
    names(sampled_gs)=paste0('null_',c(1:perm.n))
    return(sampled_gs)})
} else if (null_type='spin_brain'){
   perm.id = load_permid(data_path=sprintf('%s/BrainInfo/perm_id_weights',data_path),
                        atlas=atlas,
                        type='spin_brain')
  null_brain_data = generate_null_brain_data(brain_data, perm.id)
}

rmse.true.list=list()
rmse.null.list=list()
for (gs_idx in 1:length(geneSetList)){
gs=geneSetList[[gs_idx]]
gs_name=names(geneSetList)[gs_idx]

# true model
tmp.res=get_true_rmse(gene_data,brain_data,gs)
rmse.true=tmp.res$rmse.true
ncomp2use=tmp.res$ncomp2use

if (null_type %in% c('random_gene_coexp', 'random_gene_subset', 'random_gene')){
  sampled_gs=sampled_geneSetList[[gs_idx]]
  rmse.null=get_null_rmse_sampled_gs(gene_data, brain_data, sampled_gs, ncomp2use)
} else if (null_type=='spin_brain'){
  rmse.null=get_null_rmse_spin_brain(gene_data, null_brain_data, gs, ncomp2use)
}
rmse.true.list[[gs_name]]=rmse.true
rmse.null.list[[gs_name]]=rmse.null
}

pvals=caculate_pvals(rmse.true.list, rmse.null.list)


pvals.df=data.frame(geneSet = names(pvals), pvals.PLSR = unlist(pvals),
                    sim_brain=slrum_idx,row.names = NULL)

write.csv(pvals.df, file=savename, row.names=F)






