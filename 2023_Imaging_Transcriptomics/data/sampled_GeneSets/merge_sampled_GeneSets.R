
args = commandArgs(trailingOnly=TRUE)
data_path=args[1]
# find rds files in data_path
rds_files=list.files(data_path, pattern = "*.rds", full.names = TRUE)
# load rds files
geneSetList.final=list()
for (i in 1:length(rds_files)){
    rds_file_i=sprintf('%s/geneList_%s.rds',data_path,i) # numeric order: 1->n
  print(rds_file_i)
  tmp=readRDS(rds_file_i) 
  geneSetList.final=c(geneSetList.final,tmp)
  # remove rds_file_i
  
  print(length(geneSetList.final))
}
# get name from data_path
geneSetList_name=paste0(basename(data_path),'.rds')
saveRDS(geneSetList.final, geneSetList_name)
# copy file to parent folder
file.copy(geneSetList_name, dirname(data_path))

lapply(rds_files, file.remove)

