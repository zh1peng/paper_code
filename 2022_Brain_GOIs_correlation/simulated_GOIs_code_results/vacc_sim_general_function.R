args = commandArgs(trailingOnly=TRUE)
slrum.array.idx=as.numeric(args[1])
random_brain_type=args[2]
analysis_type=args[3]

data.path='/gpfs1/home/z/c/zcao4/Coexpression_analysis/data'
source('/gpfs1/home/z/c/zcao4/Coexpression_analysis/functions.R')

save.path='/gpfs1/home/z/c/zcao4/Coexpression_analysis/all_sim_results'
gall_expr_df=get_expression.shin2stage(data.path,hem='L')

n_random_FPR=1000

cat(sprintf(' slrum array: %s\n random brain type: %s\n analysis_type: %s\n FRP_n: %s\n',
            slrum.array.idx,random_brain_type,analysis_type,n_random_FPR))
save_name=sprintf('%s/sim_%s_%s_%s.csv',
                  save.path, random_brain_type, analysis_type,slrum.array.idx)



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




if (analysis_type=='no_coexp'){
  perm.id=NULL
} else if (analysis_type=='ctrl_coexp_nospin'){
  perm.id=readRDS(sprintf('%s/34perm_id_nospin10k.rds',data.path))
} else if  (analysis_type=='ctrl_coexp_spin') {
  perm.id=readRDS(sprintf('%s/34perm_id_spin10k.rds',data.path))
}



n.seq=seq(5, 150, 5) 
sim_seq_res=list()
for (sim_i in n.seq){
  cat(sprintf('glist_n: %s\n',sim_i))
  sim_seq_res[[sim_i]]=simulating_glist(glist_n=sim_i,
                                         sim_rep=1,
                                         gall_expr_df,
                                         SIMseed=slrum.array.idx,
                                         random_brain4FPR = random_brain_data,
                                         shuffle.id = perm.id)
}

sim_seq_res_df=do.call(rbind,sim_seq_res)
write.csv(sim_seq_res_df,save_name,row.names = F)
