args = commandArgs(trailingOnly=TRUE)
slrum.array.idx=as.numeric(args[1])
random_brain_type=args[2]


data.path='/gpfs1/home/z/c/zcao4/Coexpression_analysis/data'
source('/gpfs1/home/z/c/zcao4/Coexpression_analysis/functions.R')

save.path='/gpfs1/home/z/c/zcao4/Coexpression_analysis/SynGO_results'
gall_expr_df=get_expression.shin2stage(data.path,hem='L')

n_random_FPR=1000

cat(sprintf(' slrum array: %s\n random brain type: %s\n FRP_n: %s\n',
            slrum.array.idx,random_brain_type,n_random_FPR))
save_name=sprintf('%s/sim_%s_%s.csv',
                  save.path, random_brain_type,slrum.array.idx)



if (random_brain_type=='totally_random'){
  set.seed(1990)
  random_brain_data=matrix(runif(34*n_random_FPR,min=-1,max=1),nrow=34)
} else if (random_brain_type=='preserve_spatial'){
  set.seed(1990)
  random_brain4FPR_spin10k=readRDS(sprintf('%s/random_brain4FPR_spin10k.rds',data.path))
  idx=sample(ncol(random_brain4FPR_spin10k), n_random_FPR, replace = F)
  random_brain_data=random_brain4FPR_spin10k[,idx]
} else if (random_brain_type=='random_walk'){
    random_brain_data=readRDS(sprintf('%s/random_brain4FPR_randomwalk.rds',data.path))
} else if (random_brain_type=='RWrandom_assign300'){
    random_brain_data=readRDS(sprintf('%s/random_brain4FPR_RWrandom_assign300.rds',data.path))
} else if (random_brain_type=='uniform_subset1000'){
    random_brain_data=readRDS(sprintf('%s/random_brain4FPR_uniform_subset1000.rds',data.path))
}




perm.id.nocoexp=NULL
perm.id.nospin=readRDS(sprintf('%s/34perm_id_nospin10k.rds',data.path))
perm.id.spin=readRDS(sprintf('%s/34perm_id_spin10k.rds',data.path))



syngo_df=read.csv(sprintf('%s/syngo_int.csv',data.path))
glist=split_gene(as.character(syngo_df$int_gene[slrum.array.idx]))
glist_expr_df=get_glist_expression(gall_expr_df,glist)
glist_n=ncol(glist_expr_df)
SynGO_res=list()
SynGO_res[['GO_term']]=syngo_df$`GO term ID`[slrum.array.idx]
SynGO_res[['glist_n']]=glist_n
coexp_df=caculate_coexp(glist_expr_df)
SynGO_res[['coexp']]=coexp_df$mean
SynGO_res[['coexp_pos_adj']]=coexp_df$adj_posmean
SynGO_res[['coexp_neg_adj']]=coexp_df$adj_negmean

print('caculating FPR 1')
FPR_vect=caculate_FPR(glist_expr_df,gall_expr_df, random_brain4FPR=random_brain_data, shuffle.id=perm.id.nocoexp)
SynGO_res[['No_coexp.FPR.orig']]=FPR_vect[1]
SynGO_res[['No_coexp.FPR.sqr']]=FPR_vect[2]
SynGO_res[['No_coexp.FPR.abs']]=FPR_vect[3]
SynGO_res[['No_coexp.FPR.nsig']]=FPR_vect[4]
SynGO_res[['No_coexp.FPR.nsig.fdr']]=FPR_vect[5] 

print('caculating FPR 2')
FPR_vect=caculate_FPR(glist_expr_df,gall_expr_df, random_brain4FPR=random_brain_data, shuffle.id=perm.id.nospin)
SynGO_res[['No_spin.FPR.orig']]=FPR_vect[1]
SynGO_res[['No_spin.FPR.sqr']]=FPR_vect[2]
SynGO_res[['No_spin.FPR.abs']]=FPR_vect[3]
SynGO_res[['No_spin.FPR.nsig']]=FPR_vect[4]
SynGO_res[['No_spin.FPR.nsig.fdr']]=FPR_vect[5] 

print('caculating FPR 3')
FPR_vect=caculate_FPR(glist_expr_df,gall_expr_df, random_brain4FPR=random_brain_data, shuffle.id = perm.id.spin)
SynGO_res[['Spin.FPR.orig']]=FPR_vect[1]
SynGO_res[['Spin.FPR.sqr']]=FPR_vect[2]
SynGO_res[['Spin.FPR.abs']]=FPR_vect[3]
SynGO_res[['Spin.FPR.nsig']]=FPR_vect[4]
SynGO_res[['Spin.FPR.nsig.fdr']]=FPR_vect[5] 

SynGO_res_df=data.frame(SynGO_res)

write.csv(SynGO_res_df,save_name,row.names = F)
