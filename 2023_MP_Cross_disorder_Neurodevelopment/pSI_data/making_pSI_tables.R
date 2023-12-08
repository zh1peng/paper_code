##### make self-defined gene sets
load(sprintf('%s/pSI/pSI_SEA.Rdata',code_path))

# 
pSI=pSI %>% dplyr::select(starts_with('Cortex'))

for (pSI.threshold in c(0.1, 0.05, 0.025, 0.01,0.005,0.0025, 0.001)){
gene_list=list()
gene_n=list()

for (col2use in colnames(pSI)){
  all_gene_names=rownames(pSI)
  col_data=pSI[,col2use]
  gene_names=all_gene_names[col_data<pSI.threshold &!is.na(col_data)]
  gene_list[[col2use]]=paste(gene_names,collapse = ",")
  gene_n[[col2use]]=length(gene_names)

}
df2save=data.frame(cLabel=names(gene_list),
                   cDesc=names(gene_list),
                   cSize=unlist(gene_n),
                  cGenes=unlist(gene_list))
write.csv(df2save,sprintf('%s/pSI/pSI_only_cortex_table_%s.csv',code_path,pSI.threshold),row.names = F)
}
