# test code
code_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code'
source(file.path(code_path,'functions','data_functions.R'))
source(file.path(code_path,'functions','cor_functions.R'))
source(file.path(code_path,'functions','spinning_brain_functions.R'))
source(file.path(code_path,'functions','resampling_gene_functions.R'))

# test load_BrainDat
brain_data=load_BrainDat(atlas='desikan',type='sim_spatial0.03', col_idx = 'all')


# test load_GeneExp
gene_data=load_GeneExp(atlas='desikan',rdonor='r0.6')
b=corr_brain_gene(gene_data,brain_data,'pearson')
a=corr_brain_gene(gene_data, brain_data,'pls1w')
a.shuffled=apply(a,2,sample)
rownames(a.shuffled)=rownames(a)
hist(a[,1])

# test load_GeneSets
geneSetList=load_GeneSets(atlas = 'desikan', rdonor='r0.6',gs_type='SynGO')
geneSetList=load_GeneSets(atlas = 'desikan', rdonor='r0.4',gs_type='MF')


# test cor
geneList.cor=corr_brain_gene(gene_data, brain_data)

# test loo_pearson
geneList.corLoo=corr_brain_gene(gene_data, brain_data,method='loo')

# test pls1
geneList.pls1=corr_brain_gene(gene_data, brain_data,method='pls1')

# test resampling gene 1
geneList.null=resampling_geneList(geneList.pls1,
                    perm.n=100)


# test resampling gene 2
coexp_matrix=cor(gene_data)

# sampled_gs=sample_gs_with_coexp(coexp_matrix = coexp_matrix,
#                                 gs=geneSetList[[1]],
#                                 n_target = 10)

sampled_geneSetList=sample_geneSetList_with_constrains(geneSetList = geneSetList[1],
                                                        constrain='match_coexp',
                                                        coexp_matrix = coexp_matrix)

sampled_geneSetList1=sample_geneSetList_with_constrains(geneSetList = geneSetList[2],
                                                        constrain='match_coexp',
                                                        coexp_matrix = coexp_matrix)

combined_geneSetList <- c(sampled_geneSetList, sampled_geneSetList1)

# sampled_geneSetList=lapply(geneSetList, sample_gs_matching_coexp, coexp_matrix=coexp_matrix)

geneList.null.gs=swap_geneList(geneList.true=geneList.cor, # drop=F
                            orig_gs=geneSetList[[1]],
                            sampled_gs = sampled_geneSetList[[1]])


# toy data for swap_geneList
# set.seed(123)  # for reproducibility
# # Generate geneList.true
# num_genes <- 100  # total number of genes
# geneList.true <- matrix(1:100, ncol = 1)
# rownames(geneList.true) <- paste0("Gene", 1:num_genes)
# # Generate gs (gene set)
# num_gs <- 10  # size of gs
# gs <- sample(rownames(geneList.true), num_gs)
# # Generate sampled_gs (another gene set with the same size as gs)
# sampled_gs <- replicate(10, sample(rownames(geneList.true), num_gs), simplify = FALSE)
# result <- swap_geneList(geneList.true, gs, sampled_gs)



# test resampling within a smaller subset (more expressed in brain)
# geneList.null.gs=resampling_geneList_within_subset(geneList.true =geneList.pls1[,1,drop=F], # drop=F
#                                   gs=geneSetList[[1]],
#                                   gene_subset=sample(colnames(gene_data),500,replace = F),
#                                   perm.n=100)
gene_subset=sample(colnames(gene_data),500,replace = F)
sampled_geneSetList=sample_geneSetList_with_constrains(geneSetList = geneSetList[1:3],
                                                       constrain='gene_subset',
                                                       gene_subset = gene_subset)
geneList.null.gs=swap_geneList(geneList.true=geneList.cor, # drop=F
                               orig_gs=geneSetList[[1]],
                               sampled_gs = sampled_geneSetList[[1]])


# test spin brain
perm_id=load_permid(atlas = 'desikan',
                    type='spin_brain',
                    perm.n=5000)

null_brain_data=generate_null_brain_data(brain_data,perm_id)
geneList.null=corr_brain_gene(gene_data=gene_data,
                              brain_data=null_brain_data,
                              method='pls1')

# test aggregate gs 
gs.scores=aggregate_geneSetList(geneSetList = geneSetList,
                                geneList = geneList.pls1,
                                method='mean')

gs.scores.null=aggregate_geneSetList(geneSetList = geneSetList,
                                geneList = geneList.null,
                                method='mean')

# test aggregate gs matching coexp
gs.scores=aggregate_geneSetList_with_constrains(geneSetList = geneSetList[1:3],
                                                geneList = geneList.cor,
                                                constrain='match_coexp',
                                                coexp_matrix = coexp_matrix,
                                                method='mean')

gs.scores=aggregate_geneSetList_with_constrains(geneSetList = geneSetList[1:3],
                                                geneList = geneList.pls1,
                                                constrains='gene_subset',
                                                gene_subset = sample(colnames(gene_data),500),
                                                method='mean')


###########################
brain_data=load_BrainDat(csvname='desikan_sim_spatial0.03',col_idx = 1)
gene_data=load_GeneExp(atlas='desikan',rdonor='r0.6')
geneSetList=load_GeneSets(atlas = 'desikan', rdonor='r0.6',gs_type='SynGO')
geneList.true=corr_brain_gene(gene_data, brain_data)
perm_id=load_permid(atlas = 'desikan',
                    type='spin_brain',
                    perm.n=5000)

null_brain_data=generate_null_brain_data(brain_data,perm_id)
geneList.null=corr_brain_gene(gene_data=gene_data,
                              brain_data=null_brain_data)
cat('==========\n')
cat('start testing \n')
tt=Sys.time()
res=list()
#for (method2test in c('mean','median')){
  cat(paste0(method2test,'....\n'))
  score.true=aggregate_geneSetList(geneSetList, geneList.true,method=method2test)
  score.null=aggregate_geneSetList(geneSetList[1:3], geneList.null,method='mean')
  
  score.null=aggregate_geneSetList_with_constrains(geneSetList = geneSetList[1:3],
                                                    geneList = geneList.true,
                                                    constrain='gene_subset',
                                                    gene_subset = sample(colnames(gene_data),500),
                                                    method='mean')
  
  pvals=caculate_pvals(score.true,score.null,method='standard')
#  res[[method2test]][['pvals']]=pvals
#}


df.list=list()
for (var2extract in c('pvals')){
  df.list[[var2extract]]= as.data.frame(lapply(res, function(x){do.call(rbind,x[[var2extract]])}))
}
df2save=do.call(cbind,df.list) %>%  tibble::rownames_to_column("geneSet") %>% mutate(sim_brain=1)

Sys.time()-tt


