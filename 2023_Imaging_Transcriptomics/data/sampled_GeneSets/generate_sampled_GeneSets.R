args = commandArgs(trailingOnly=TRUE)

atlas=args[1] # desikan, schaefer100, schaefer200
rdonor=args[2] # r0.2, r0.4, r0.6
gs_type=args[3] # MF, Sim, SynGO
constrain=args[4] # match_coexp, gene_subset
slrum_idx=as.numeric(args[5]) # 1,2,3,4,5

#print input
cat('==========\n')
cat(sprintf('atlas: %s\n',atlas))
cat(sprintf('rdonor: %s\n',rdonor))
cat(sprintf('gs_type: %s\n',gs_type))
cat(sprintf('constrain: %s\n',constrain))
cat(sprintf('slrum_idx: %s\n',slrum_idx))
cat('==========\n')

library(dplyr)

data_path='/gpfs1/home/z/c/zcao4/revision_code/data'


code_path='/gpfs1/home/z/c/zcao4/revision_code/functions'

# source functions
source(sprintf( '%s/data_functions.R',code_path))
source(sprintf( '%s/resampling_gene_functions.R',code_path))


geneSetList=load_GeneSets(data_path=sprintf('%s/GeneSets',data_path),
                        atlas=atlas,
                        rdonor=rdonor,
                        gs_type=gs_type)
gene_data=load_GeneExp(data_path=sprintf('%s/GeneExp',data_path),
                        atlas=atlas,rdonor=rdonor)
# matching coexp
if (constrain=='match_coexp'){

save_path=sprintf('%s/sampled_GeneSets/%s_%s_%s_%s',data_path,atlas,rdonor,gs_type, constrain)
#create save path
if (!dir.exists(save_path)){
  dir.create(save_path)
}
savename=sprintf('%s/geneList_%s.rds',save_path, slrum_idx)

coexp_matrix=cor(gene_data)

gene_set_size=length(geneSetList[[slrum_idx]])
gene_set_name=names(geneSetList)[slrum_idx]

tt=Sys.time()
tmp_geneSetList=sample_geneSetList_with_constrains(geneSetList = geneSetList[slrum_idx],
                                                        constrain='match_coexp',
                                                        coexp_matrix = coexp_matrix)
ttt=Sys.time()-tt

gs_coexp_matrix <- coexp_matrix[geneSetList[[slrum_idx]], geneSetList[[slrum_idx]]]
gs_coexp_lower <- gs_coexp_matrix[lower.tri(gs_coexp_matrix)]
target_coexp <- mean(gs_coexp_lower)

all_coexp=mean(coexp_matrix[lower.tri(coexp_matrix)])

saveRDS(tmp_geneSetList,savename)

# cat(sprintf('time: %s',ttt))
# cat(sprintf('gene set name: %s',gene_set_name))
# cat('gene set size: %s', gene_set_size)

# append gene_set_name, gene_set_size, time to a file
#append_file=sprintf('%s/task_info.csv',save_path)
#write.csv(data.frame(gene_set_name, gene_set_size, ttt), append_file, row.names = F, append = T)

data = paste(gene_set_name, gene_set_size, target_coexp, ttt, all_coexp, sep = ",")
write(data, file = sprintf('%s/task_info.csv', save_path), append = TRUE)

} else if(constrain=='gene_subset'){
  brain_body_info=read.csv(sprintf('%s/sampled_GeneSets/Background_genes_null_brain.csv',data_path))%>% filter(`brain...body..default.`==1)
  gene_subset=base::intersect(brain_body_info$symbols,colnames(gene_data))
  sampled_geneSetList=sample_geneSetList_with_constrains(geneSetList = geneSetList,
                                                constrain='gene_subset',
                                                gene_subset = gene_subset)
  saveRDS(sampled_geneSetList,sprintf('%s/sampled_GeneSets/%s_%s_%s_%s.rds',data_path, atlas,rdonor,gs_type, constrain))

}


