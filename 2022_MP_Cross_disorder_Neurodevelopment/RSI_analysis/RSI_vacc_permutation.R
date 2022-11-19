library(dplyr)
library(tidyr)

source('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/code/pca_functions.R')
source('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/code/external_functions.R')


# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
slurm_arrayid <- as.numeric(slurm_arrayid)
loop_start=(slurm_arrayid-1)*100+1
loop_end=slurm_arrayid*100

time_name='T2'
rsi_name='nd'
rsi_df=read.csv(sprintf('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/ABCD_RSI/ABCD_%s_%s_ready.csv',time_name, rsi_name))
rsi_df=preprocess_df(rsi_df)


# load CT
ct_df=read.csv(sprintf('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/ABCD_CT/ABCD_%s_CT_ready.csv',time_name)) %>% filter(SubjID %in% rsi_df$SubjID)
ct_df=preprocess_df(ct_df)
ct.pca=do_pca_df(ct_df,if_resid = F)
ct.pca=keep_pc1(ct.pca)

pc1_mat=t(as.matrix(ct.pca$pc1_only_data))
#rsi_mat=t(as.matrix(rsi_df.resid %>% select(starts_with('resid_L_')|starts_with('resid_R_'))))
#rsi_mat=t(as.matrix(rsi_df %>% select(starts_with('L_')|starts_with('R_')))) # no idea why This is not woking on vacc
rsi_mat=t(as.matrix(rsi_df[,colnames(rsi_df)[grep('^L_|^R_',colnames(rsi_df))]]))


# shuffle PC1 data for each participants [ need to spin to create null]
perm.id=readRDS('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/code/68perm_id.rds')
all_null_r=numeric()

for (perm_i in c(loop_start:loop_end)){
  print(perm_i)
  shuffle.idx=perm.id[,perm_i]
  pc1_mat.shuffled=pc1_mat[shuffle.idx,]
  all_null_r=c(all_null_r,mean(diag(cor(pc1_mat.shuffled,rsi_mat,method = 'pearson'))))
}
write.table(all_null_r, file=sprintf('/gpfs1/home/z/c/zcao4/PCA_project/RSI_permutation/Null_dist_PC1_RSI/%s_%s_null.csv',
                                     time_name,rsi_name),
            append = T,row.names = F,col.names = F)

