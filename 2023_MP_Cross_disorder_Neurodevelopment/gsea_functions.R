# import functions from DOSE and clusterProfiler
gseaScores=getFromNamespace('gseaScores','DOSE')
geneSet_filter=getFromNamespace('geneSet_filter','DOSE')
get_geneSet_index=getFromNamespace('get_geneSet_index','DOSE')
getGeneSet=getFromNamespace('getGeneSet','DOSE')
get_GO2TERM_table=getFromNamespace('get_GO2TERM_table','clusterProfiler')
get_GOTERM=getFromNamespace('get_GOTERM','clusterProfiler')
build_Anno=getFromNamespace("build_Anno", "DOSE")
get_GO_Env=getFromNamespace('get_GO_Env','clusterProfiler')
get_GO_data=getFromNamespace('get_GO_data','clusterProfiler')
build_Anno <- getFromNamespace("build_Anno", "DOSE")
calculate_qvalue <- getFromNamespace('calculate_qvalue','DOSE')
TERM2NAME <- getFromNamespace('TERM2NAME','DOSE')
GSEA_internal <- getFromNamespace('GSEA_internal','clusterProfiler')


# map SYMBOL to ENTREZID
convert_gene_ID <- function(r_vals){
  # map SYMBOL to ENTREZID
  # convert r_vals to r_df
  r_df=as.data.frame(r_vals) %>% 
    mutate(SYMBOL=sub('gene_','',rownames(.))) %>% 
    # will use SYMBOL for the analysis/keep consistant with GCSA
    # mutate(ENTREZID=mapIds(x=org.Hs.eg.db,
    #                        keys=SYMBOL,
    #                        keytype = "SYMBOL",
    #                        column = "ENTREZID")) %>%
    # filter(complete.cases(ENTREZID))
  return(r_df)
}

# make a geneList for the subsequent analysis
make_geneList <- function(r_df,
                          colname2use='brain_PC1'){
  # make a geneList for the subsequent analysis
  # based on the colname
  # brain_PC1 is the true correlation
  # null_brain_* is the null correlation
  geneList=r_df[,colname2use]
  names(geneList)=as.character(r_df[,'SYMBOL']) # will use SYMBOL for all analysis
  geneList = sort(geneList, decreasing = TRUE)
  return(geneList)
}

# caculate the enrichment score  
get_es <- function(selected.gs,geneList){
  # caculate the enrichment score for the selected.gs and input geneList
  observed_info <- lapply(selected.gs, function(gs)
    gseaScores(geneSet=gs,
               geneList=geneList,
               exponent=1))
  observedScore <- sapply(observed_info, function(x) x$ES)
  return(observedScore)
}


#library(doParallel)
# function for parallel function
# get_null_es_single <- function(null_r_mx_i,selected.gs){
#       null_geneList=sort(null_r_mx_i,decreasing = TRUE)
#       es=get_es(selected.gs,null_geneList)
#       return(es)
#       }


get_es_batch <- function(r_vals,
                         selected.gs,
                         if_r2z=T, # convert r to z vals Fisher
                         es_type=c('true_model','null_model'), 
                         null_iter=10000){
  if (if_r2z==T){
    r_vals=apply(r_vals,2,function(r){0.5*(log((1+r)/(1-r)))})
  }
  
  if (es_type=='true_model'){
    true_r_df=convert_gene_ID(r_vals)
    true_geneList=make_geneList(true_r_df,'brain_PC1')
    true_es=get_es(selected.gs,true_geneList)
    es_df=data.frame(true_es)
  } else if (es_type=='null_model'){
    
    # FOR LOOP for 5 null model: 29s
    null_r_df=convert_gene_ID(r_vals) 
    pb = txtProgressBar(min = 0, max = null_iter, initial = 0, style = 3)
    null_es_list=list() 
    for (null_i in c(1:null_iter)){
       null_colname=paste0('null_brain_',null_i)
       null_geneList=make_geneList(null_r_df,null_colname)
       null_es_list[[null_colname]]=get_es(selected.gs, null_geneList)
       setTxtProgressBar(pb,null_i)
     }
    es_df=data.frame(null_es_list)
  }
    
    ### sapply for 5 null model: 29s
    # null_r_df=convert_gene_ID(r_vals) 
    # rownames(null_r_df)=as.character(null_r_df[,'ENTREZID']) # set rownames
    # null_r_df=null_r_df%>% dplyr::select(starts_with('null_brain_')) # select data
    # null_r_mx=as.matrix(null_r_df)
    # ptm <- proc.time()
    # null_es=sapply(c(1:5),function(x,...){
    #   print(x)
    #   null_geneList=sort(null_r_mx[,x],decreasing = TRUE)
    #   es=get_es(selected.gs,null_geneList)
    # })
    #
  return(es_df)
  }
  




get_pvals <- function(true_df,null_df){
  print(sprintf('calculating pvals for %d GO terms', dim(true_df)[1]))
  pb = txtProgressBar(min = 0, max = dim(true_df)[1], initial = 0, style = 3)
  pvals=matrix(NA,nrow=dim(true_df)[1],ncol=1)
  for (row_i in c(1:dim(true_df)[1])){
    if( is.na(true_df[row_i,1]) ) {
      pvals[row_i]=NA
    } else if (true_df[row_i,1]>= 0) {
      pvals[row_i,1]=(sum(null_df[row_i, ] >= true_df[row_i,1]) +1) / (sum(null_df[row_i,] >= 0) +1)
    } else { # true_es_df[i,] < 0
      pvals[row_i,1]=(sum(null_df[row_i, ] <= true_df[row_i,1]) +1) / (sum(null_df[row_i,] < 0) +1)
    }
    setTxtProgressBar(pb,row_i)
  }
  rownames(pvals)=rownames(true_df)
  return(pvals)
}

# a slightly different way to calculate p values, give slightly different results.
# get_pvals2 <- function(true_df,null_df){
#   print(sprintf('calculating pvals for %d GO terms', dim(true_df)[1]))
#   pb = txtProgressBar(min = 0, max = dim(true_df)[1], initial = 0, style = 3)
#   pvals=matrix(NA,nrow=dim(true_df)[1],ncol=1)
#   for (row_i in c(1:dim(true_df)[1])){
#     mean.under.null=mean(null_df[row_i,])
#     null.aug=c(true_df[row_i,1],null_df[row_i,])
#     null.aug.adj=null.aug-mean.under.null
#     true.adj=true_df[row_i,1]-mean.under.null
#     pvals[row_i,]=sum(abs(null.aug.adj) >=abs(true.adj))/length(null.aug.adj)
#     setTxtProgressBar(pb,row_i)
#   }
#   rownames(pvals)=rownames(true_df)
#   return(pvals)
# }
  
  



normalized_ES <- function(ES, pos.m, neg.m) {
  s <- sign(ES)
  m <- numeric(length(ES))
  m[s==1] <- pos.m[s==1]
  m[s==-1] <- neg.m[s==-1]
  ES/m
}


# GCEA method
gceaScores <- function(geneSet,
                       df_val,
                       method2use=c('mean','median','absmean'),
                       ID2use='SYMBOL'){
  geneSet = intersect(geneSet, df_val[,ID2use])
  hits = df_val[,ID2use] %in% geneSet ## logical
  df_val.hit= df_val[hits,] %>% dplyr::select(-c('SYMBOL'))
  
  if (method2use=='mean'){
    CS=apply(df_val.hit, 2, mean)  
  } 
  else if (method2use=='median'){
    CS=apply(df_val.hit, 2, median)  
  }
  else if (method2use=='absmean'){
    CS=apply(df_val.hit,2,function(x){mean(abs(x))})
  }
  return(CS)
}

# caculate the enrichment score  
get_cs <- function(r_vals,
                   selected.gs,
                   if_r2z=T,
                   method2use='mean',
                   ID2use='SYMBOL'){
  if (if_r2z==T){
    r_vals=apply(r_vals,2,function(r){0.5*(log((1+r)/(1-r)))})
  }
  r_df=convert_gene_ID(r_vals) # this creates r_df; SYMBOL was used
  # caculate the enrichment score for the selected.gs and input geneList
  observed_info <- lapply(selected.gs, function(gs)
    gceaScores(geneSet=gs,
               df_val=r_df,
               method2use=method2use,
               ID2use=ID2use))
  observedScore = t(data.frame(observed_info))
  rownames(observedScore)=names(observed_info) # this is to fix ':' became '.' in ID
  return(observedScore)
}



create_USER_DATA <- function(GS_file){
GS_df=read.csv(GS_file, stringsAsFactors=FALSE)
TERM2GENE=GS_df%>%dplyr::select(cLabel, cGenes) %>% # TERM2GENE
            dplyr::rename(geneID=cGenes) %>% 
            dplyr::mutate(geneID = strsplit(geneID, ',')) %>%
            tidyr::unnest(cols = c(geneID))  
TERM2NAME=GS_df%>%dplyr::select(cLabel, cDesc) # cDesc
USER_DATA = build_Anno(TERM2GENE,TERM2NAME)
 return(USER_DATA)     
}

# EXAMPLE: 
# GS_file='F:\\Google Drive\\post-doc\\p_factor\\ABAnnotate-main\\datasets\\DisGeNET-diseaseCuratedMental-discrete.csv'
# USER_DATA=create_USER_DATA(GS_file)

# USER_DATA can be loaded using org.Hs.eg.db as well
# go_env=get_GO_data('org.Hs.eg.db',ont2use,'ENTREZID') # or 'SYMBOL'
# USER_DATA=getGeneSet(go_env)



# a general function to do enrichment analysis
# USER_DATA is the geneSets structure used for clusterProfiler
enrich_with_snull <- function(true_r_vals,
                              null_r_vals,
                              USER_DATA,
                              minGSSize=10,
                              maxGSSize=200,
                              by=c('gsea','gcea'),
                              returnDF=T){
  if(missing(by)){by='gsea'}
# create true_geneLIst for: 
# 1. filter GS size. 
# 2. get the information using gsea_internal 
true_z_vals=apply(true_r_vals,2,function(r){0.5*(log((1+r)/(1-r)))})
true_z_df=convert_gene_ID(true_z_vals)
true_geneList=make_geneList(true_z_df,'brain_PC1')

# filter GS size [minGSSize,maxGSSize]
USER_GS=getGeneSet(USER_DATA)
# minGSSize=10
# maxGSSize=200
selected.gs=geneSet_filter(USER_GS, true_geneList, minGSSize,maxGSSize)

if (by=='gsea'){
#################################
# caculate ES for true_r_vals
true_es_df=get_es_batch(true_r_vals,
                        selected.gs,
                        if_r2z=T, 
                        es_type='true_model')

# caculate ES for null_r_vals
null_es_df=get_es_batch(null_r_vals,
                        selected.gs,
                        if_r2z=T, 
                        es_type='null_model',
                        null_iter = dim(null_r_vals)[2]) #dim(null_r_vals)[2]
# calculate NES
pos.m <- apply(null_es_df, 1, function(x) mean(x[x >= 0]))
neg.m <- apply(null_es_df, 1, function(x) abs(mean(x[x < 0])))
true_df <- normalized_ES(true_es_df, pos.m, neg.m)
null_df <- apply(null_es_df, 2, normalized_ES, pos.m=pos.m, neg.m=neg.m)
} else if (by=='gcea'){
#################################
# caculate CS for true_r_vals
  true_df=get_cs(true_r_vals,
                 selected.gs,
                 if_r2z=T,
                 method2use='mean',
                 ID2use = 'SYMBOL')
  
# caculate CS for null_r_vals
  null_df=get_cs(null_r_vals,
                 selected.gs,
                 if_r2z=T,
                 method2use='mean',
                 ID2use = 'SYMBOL')
  # no need to adjust for size
  # pos.m <- apply(null_df, 1, function(x) mean(x[x >= 0]))
  # neg.m <- apply(null_df, 1, function(x) abs(mean(x[x < 0])))
  # true_df <- normalized_ES(true_df, pos.m, neg.m)
  # null_df <- apply(null_df, 2, normalized_ES, pos.m=pos.m, neg.m=neg.m)
  
}

#################################
# get pvals
# make sure the rows are in the same order
if (!identical(rownames(true_df),rownames(null_df))){
  stop("rows in true_df and null_df are NOT in the same order")
}

p_df=get_pvals(true_df, null_df)
pvals.adj=p.adjust(p_df,'BH')# fdr adjusted  pvals
qvals = calculate_qvalue(p_df)# q values

# save pvals information
if (by=='gsea'){
sp_null_df=cbind(true_es_df,true_df,p_df,pvals.adj,qvals)
colnames(sp_null_df)=c('sn_ES','sn_NES','sn_NESpval','sn_NESp.adjust','sn_NESqvals')
}
else if (by=='gcea'){
sp_null_df=cbind(true_df,p_df,pvals.adj,qvals)
colnames(sp_null_df)=c('sn_CS','sn_CSpval','sn_CSp.adjust','sn_CSqvals')
}

# get other basic information
res=GSEA_internal(geneList    = true_geneList,
                exponent      = 1,# weighted by the correlation value
                minGSSize     = minGSSize,
                maxGSSize     = maxGSSize,
                eps           = 1e-10, 
                pvalueCutoff  = 1,
                pAdjustMethod = 'BH',
                USER_DATA     = USER_DATA, 
                seed          = 1990,
                verbose       = TRUE,
                by            = 'fgsea')
res_df.orig=as.data.frame(res)
res_df.merged=merge(res_df.orig,sp_null_df,by='row.names')
if(returnDF){
  return(res_df.merged)}
else if (!returnDF){
  returnlist=list()
  returnlist[['enrich_obj']]=res
  returnlist[['merged_res_df']]=res_df.merged
  return(returnlist)}
}


# function to update the returnlist for gcea
update_gcea_res.obj <- function(gcea_res_list,
                                settype='BP',
                                filter_fdrp=T, #otherwise use pvalue
                                threshold=0.05
                                ){
  updated_res=gcea_res_list$merged_res_df %>% 
    dplyr::select(-c('enrichmentScore', 'NES', 'pvalue','p.adjust',
                     'qvalues','Row.names','core_enrichment',
                     'leading_edge','rank')) %>% 
    dplyr::rename(CategoryScore=sn_CS,
                  pvalue=sn_CSpval,
                  p.adjust=sn_CSp.adjust,
                  qvalues=sn_CSqvals) 
  if (filter_fdrp){
    updated_res=updated_res[updated_res$p.adjust<threshold,]
  } else {
    updated_res=updated_res[updated_res$pvalue<threshold,]
  }

  updated_res=updated_res[order(updated_res$pvalue),]
  
  # in tree plot, the size of the node is scaled by the size of core_enrichment
  # Here create the core_enrichment colume that reflect the rank of the pvalue
  updated_res$core_enrichment=sapply(c(nrow(updated_res):1), function(x){paste0(rep('genex',x),collapse='/')})
  #updated_res$core_enrichment='genex'
  
  res.obj=gcea_res_list$enrich_obj
  res.obj@result=updated_res
  res.obj@organism='Homo sapiens'
  res.obj@setType=settype
  res.obj@keytype='SYMBOL'
  return(res.obj)
}








