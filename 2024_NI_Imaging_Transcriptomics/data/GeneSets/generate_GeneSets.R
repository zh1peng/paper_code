
# add load functions
source('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/functions/data_functions.R')

library(tidyr)
# load helper functions from clusterProfiler
get_GO_data=getFromNamespace('get_GO_data','clusterProfiler')
geneSet_filter=getFromNamespace('geneSet_filter','DOSE')
build_Anno=getFromNamespace("build_Anno", "DOSE")
getGeneSet=getFromNamespace('getGeneSet','DOSE')

simulate_geneSetList <- function(gene_data,
                                 sim.size = seq(20, 200, by=20),
                                 sim.rep = 50) {
  
  # create simulated geneSetList for each [sim.size]
  # and this is repeated for [sim.rep] times with different seeds
  lapply(sim.size,function(size_i){ 
  set.seed(size_i)
  seeds = sample(100000:999999, size = sim.rep, replace = F) # seeds for reps
  names(seeds)=paste0('S',seeds)
  lapply(seeds, function(rep_i) { # alternative: sapply: simplify=F, use.names=T
    set.seed(rep_i) # set seeds for reps
    sample(colnames(gene_data), size_i, replace = F)
  })}) %>% purrr::flatten()
  }





get_GO_geneSetList <- function(gene_data,
                               ont2use = c('MF', 'BP', 'CC'),
                               minGSSize = 20,
                               maxGSSize = 200) {

GO_env=get_GO_data('org.Hs.eg.db',ont=ont2use,'SYMBOL')
USER_DATA=getGeneSet(GO_env)

# create a tmp geneList for filtering
brain_data=rnorm(nrow(gene_data)) # create a fake brain data
geneList.tmp=cor(gene_data, brain_data)
geneSetList=geneSet_filter(USER_DATA, geneList.tmp[,1], minGSSize=minGSSize,maxGSSize=maxGSSize)
return(geneSetList)
}

get_synGO_geneSetList <- function(gene_data,
                                  synGO_file='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/syngo_ontologies.xlsx',
                                  minGSSize = 20,
                                  maxGSSize = 200) {
#check if synGO_file exist
if (!file.exists(synGO_file)) {
  stop(sprintf('synGO_file %s does not exist',synGO_file))
}

GS_df=readxl::read_xlsx(synGO_file)
TERM2GENE=GS_df%>%dplyr::select(`GO term ID`, `genes - hgnc_symbol`) %>% # TERM2GENE
                  dplyr::rename(cLabel=`GO term ID`,
                                geneID=`genes - hgnc_symbol`) %>% 
                  mutate(geneID = strsplit(geneID, ';')) %>%
                  unnest(cols = c(geneID))  
TERM2NAME=GS_df%>%dplyr::select(`GO term ID`, `GO term name`) %>% 
                  dplyr::rename(cLabel=`GO term ID`,
                                description=`GO term name`)
USER_DATA = build_Anno(TERM2GENE,TERM2NAME)
USER_GS=getGeneSet(USER_DATA)
brain_data=rnorm(nrow(gene_data)) # create a fake brain data
geneList.tmp=cor(gene_data, brain_data)
geneSetList=geneSet_filter(USER_GS, geneList.tmp[,1], minGSSize=minGSSize,maxGSSize=maxGSSize)
return(geneSetList)
}


# create geneSetList
for (atlas_i in c('desikan','schaefer100','schaefer200')){
  for (rdonor in c('r0.2','r0.4','r0.6')){
    gene_data=load_GeneExp(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneExp',
                           atlas=atlas_i,rdonor=rdonor)

    # Random500
    geneSetList=simulate_geneSetList(gene_data)
    saveRDS(geneSetList,sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/%s_%s_Sim.rds',atlas_i,rdonor))
    
    # GO-MF
    geneSetList=get_GO_geneSetList(gene_data,ont2use='MF',minGSSize=20,maxGSSize=200)
    saveRDS(geneSetList,sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/%s_%s_MF.rds',atlas_i,rdonor))
    
    # SynGO
    geneSetList=get_synGO_geneSetList(gene_data,
                                      synGO_file='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/syngo_ontologies.xlsx',
                                      minGSSize=20,maxGSSize=200)
    saveRDS(geneSetList,sprintf('F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data/GeneSets/%s_%s_SynGO.rds',atlas_i,rdonor))
  }
}




