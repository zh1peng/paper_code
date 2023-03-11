merge_sim_results<-function(sim_path,
                            sim_n=1000,
                            remove.orig=T){
  files.names<-list.files(sim_path,pattern="sim_.*\\.csv",full.names = T) # reg exp
  print(files.names)
  if(length(files.names)==sim_n){ # do following if the length of file meet a certain number
  # read csv files and return them as items in a list()
  theList = lapply(files.names, read.csv, stringsAsFactors=F)
  theResult = do.call(rbind,theList)
  write.csv(theResult, sprintf('%s/Res_%s_sim%s.csv',sim_path,basename(sim_path),length(theList)), row.names = F)
  if (remove.orig){
    tmp=lapply(files.names, file.remove)
  }
  }
}
#save_path='F:/Google Drive/post-doc/vitural_histology_revisit/new_code/real_res'
save_path='/gpfs1/home/z/c/zcao4/geneset_sim_analysis/all_sim'
all_path=list.dirs(save_path,recursive = F)
lapply(all_path,merge_sim_results,sim_n=7)

# this is by geneSet
get_psig <- function(res.file,
                     threshold = 0.05) {
  # caculate the psig across random maps
  
  # read
  res.df = read.csv(res.file, stringsAsFactors = F)
  
  # select start with pvals and create adj.pvals
  pvals.df=res.df %>%
    select(c(starts_with('pvals'), geneSet))

  pvals.nested = pvals.df %>%
    group_by(geneSet) %>%
    nest() %>%
    mutate(pvals.fdr = purrr::map(data, function(x) {apply(x, 2, p.adjust, method = 'fdr')}))

  psig.df = pvals.nested %>%
    mutate(psig = purrr::map(data, function(x) {
      colMeans(x < 0.05)
    })) %>%
    select(-data, -pvals.fdr) %>%
    unnest_wider(psig) %>% # use unnest_wider
    ungroup() %>%
    rename_all( ~ stringr::str_replace(., '^pvals.', ''))
  
  psig.df_fdr = pvals.nested %>%
    mutate(psig = purrr::map(pvals.fdr, function(x) {colMeans(x < 0.05)})) %>%
    select(-data, -pvals.fdr) %>%
    unnest_wider(psig) %>% # use unnest_wider
    ungroup() %>%
    rename_all( ~ stringr::str_replace(., '^pvals.', ''))

  return(list(pvals=pvals.df,
              psig=psig.df,
              psig_fdr=psig.df_fdr))
}

extract_psig_df <- function(sim_path){
  res.file<-list.files(sim_path,pattern="*.csv",full.names = T)
  if (length(res.file)==1){
  tmp=get_psig(res.file[1], threshold = 0.05)
  write.csv(tmp$pvals,sprintf('%s/Pvals_%s.csv',sim_path,basename(sim_path)), row.names = F)
  write.csv(tmp$psig,sprintf('%s/Psig_%s.csv',sim_path,basename(sim_path)), row.names = F)
  write.csv(tmp$psig_fdr,sprintf('%s/FDR_Psig_%s.csv',sim_path,basename(sim_path)), row.names = F)
  }
}

save_path='/gpfs1/home/z/c/zcao4/geneset_sim_analysis/all_sim'
all_path=list.dirs(save_path,recursive = F)
lapply(all_path,extract_psig_df)



  




