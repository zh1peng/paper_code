library(dplyr)
# read gene expression data
get_expression.shin2stage <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data',
                                      hem='L'){
  dk_labels=read.csv(sprintf('%s/dk_standard_var_order.csv',data_path),stringsAsFactors = F)
  dk_labels.hem=dk_labels %>% filter(grepl(sprintf('^%s_',hem), label))
  # read gene expression data
  gall_expr_df = read.table(sprintf("%s/AllenHBA_DK_ExpressionMatrix.tsv", data_path),
                            stringsAsFactors = F)
  gall_expr_df=as.data.frame(t(gall_expr_df))
  rownames(gall_expr_df) = gsub("ctx.lh.", "L_", rownames(gall_expr_df))
  rownames(gall_expr_df) = gsub("ctx.rh.", "R_", rownames(gall_expr_df))
  ## Load gene symbols and other information for the reference panel 
  filtered_gene_df = read.csv(sprintf('%s/Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv',data_path),
                              sep = '\t',stringsAsFactors = F)
  gall_expr_df=gall_expr_df %>% filter(rownames(.) %in% dk_labels.hem$label) %>% 
    select(filtered_gene_df$GeneSymbol)
  match.idx=match(dk_labels.hem$label,rownames(gall_expr_df))
  gall_expr_df=gall_expr_df[match.idx,]
  
  # file name: Shin2stage/Abgen
  #rename and check the order
  return (gall_expr_df)
}

# read brain data
read_brain_phenotype <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data',
                                 file_name='sample_brain_phenotype.csv',
                                 hem='L'){
  dk_labels=read.csv(sprintf('%s/dk_standard_var_order.csv',data_path),stringsAsFactors = F)
  dk_labels.hem=dk_labels %>% filter(grepl(sprintf('^%s_',hem), label))
  # rename and check the order
  brain_df=read.csv(sprintf('%s/%s',data_path,file_name))
  brain_df=brain_df %>% tibble::column_to_rownames(var = "label") %>% 
    filter(rownames(.) %in% dk_labels.hem$label)
  colnames(brain_df)='x'
  match.idx=match(dk_labels.hem$label,rownames(brain_df))
  brain_df=brain_df[match.idx,,drop = FALSE] #drop to one-dimension
  return (brain_df)
}



# get gene expression for a list of genes
get_glist_expression <- function(gall_expr_df, 
                                 glist){
  glist.int=intersect(glist, colnames(gall_expr_df))
  glist_expr_df=gall_expr_df %>% select(all_of(glist.int))
  return (glist_expr_df)
}




caculate_coexp <- function (glist_expr_df,
                            stattype='true'){
  coexp.matrix=cor(glist_expr_df,method = 'pearson')
  coexp.matrix.upper=coexp.matrix[upper.tri(coexp.matrix, diag = FALSE)]
  # caculate abs mean/mean/pos mean/neg mean
  # ngenes=length(genes.name)
  # ncor=(ngenes*(ngenes-1))/2
  ncor=length(coexp.matrix.upper)
  coexp.list=list()
  coexp.list[['stattype']]=stattype
  coexp.list[['ngenes']]=dim(glist_expr_df)[2]
  coexp.list[['mean']]= mean(coexp.matrix.upper) # average coexp
  coexp.list[['absmean']]= mean(abs(coexp.matrix.upper))# average abs coexp
  coexp.list[['posmean']]= mean(coexp.matrix.upper[coexp.matrix.upper>0])# average pos coexp
  coexp.list[['negmean']]= mean(coexp.matrix.upper[coexp.matrix.upper<0])# average neg coexp
  coexp.list[['possum']]= sum(coexp.matrix.upper[coexp.matrix.upper>0])# sum of pos coexp
  coexp.list[['negsum']]= sum(coexp.matrix.upper[coexp.matrix.upper<0])# sum of neg coexp
  coexp.list[['adj_posmean']]= coexp.list[['possum']]/ncor # sum of pos coexp/total number of coexp
  coexp.list[['adj_negmean']]= coexp.list[['negsum']]/ncor# sum of neg coexp/total number of coexp
  coexp_df=data.frame(coexp.list)
  return(coexp_df)
}



# Moran Idx
# data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data'
# lh_centroid_data=read.csv(sprintf('%s/lh_centroid.csv',data_path))
# centroid.l=lh_centroid_data %>%select(starts_with('centroid'))
# rownames(centroid.l)=lh_centroid_data$Row
# dist.mat=as.matrix(dist(centroid.l))
# w=1/dist.mat
# diag(w)= 0
# saveRDS(w, file = "F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data/Moran_w.rds")


get_df_moran_idx <- function(df){
  moran_idx=apply(df, 2, get_moran_idx, weight=.GlobalEnv$w) # read w first
  return(as.data.frame(moran_idx))
}

get_moran_idx <- function(x,weight){
  temp=Moran.I(x,weight)
  return(temp$observed)
}




build_null_brain_df <- function(perm.id,brain_df){
  # before put perm.id in the function 
  # need to make sure no duplicate in perm.id.
  brain_data=as.matrix(brain_df)
  perm.n=dim(perm.id)[2]
  all_null_brain=list()
  for (idx in c(1:perm.n)){
    tmp_brain_data=brain_data
    idx_name=sprintf('null_brain_%d',idx)
    perm.idx=perm.id[,idx]
    null_brain=tmp_brain_data[perm.idx]
    all_null_brain[[idx_name]]=null_brain
  }
  null_brain_df=as.data.frame(all_null_brain)
  return(null_brain_df)
}

cor2p <- function(r,n=34){
  #n is 34 here but this is how to find pairwise n
  #n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
  # a matrix version of cor.test https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}

sigtest<- function(brain_df,
                  glist_expr_df,
                  gall_expr_df,
                  shuffle.id=NULL){
  # shuffle.id determines how to shuffle the brain data (spin or nospin)
  cor.true=cor(glist_expr_df,brain_df,method='pearson')
  
  if (is.null(shuffle.id)){ # resampling approach
    
    cor.all=cor(gall_expr_df,brain_df,method='pearson')
    # re-sample r values for  random genes with same size
    size.glist=dim(glist_expr_df)[2]
    size.gall=dim(gall_expr_df)[2]
    cor.null= matrix(NA, nrow=size.glist, ncol=10000)
    for(nreps_i in 1:10000){
      row.ind = sample(size.gall,size = size.glist,replace = F)
      cor.null[,nreps_i] <- cor.all[row.ind] # cor.all[,row.ind] as it is a 1*n vector 
    } 
  } else {
  null_brain_df=build_null_brain_df(shuffle.id, brain_df)
  cor.null=cor(glist_expr_df,null_brain_df,method='pearson')  
  }

  
  # test significant number: p<0.05 or p.fdr<0.05
  p.true=cor2p(cor.true)
  # cor.test(glist_expr_df[,1],brain_df[,1])
  # p.true[1]
  p.null=cor2p(cor.null)
  # cor.test(glist_expr_df[,1],null_brain_df[,1])
  # p.null[1,1]
  # cor.test(glist_expr_df[,2],null_brain_df[,1])
  # p.null[2,1]
  #p.adjust(cor2p(cor(glist_expr_df,null_brain_df[,940])),'fdr')
  
  p.true.fdr=apply(p.true,2,p.adjust,method='fdr')
  p.null.fdr=apply(p.null,2, p.adjust,method='fdr')
  
  
  nsig.true=apply(p.true<0.05, 2, sum)
  nsig.null=apply(p.null<0.05,2,sum)
  prctile.nsig=quantile(nsig.null, c(0.95))
  if_sig.nsig=nsig.true>prctile.nsig[1]
  
  
  nsig.true.fdr=apply(p.true.fdr<0.05, 2, sum)
  nsig.null.fdr=apply(p.null.fdr<0.05, 2, sum)
  prctile.nsig.fdr=quantile(nsig.null.fdr, c(0.95))
  if_sig.nsig.fdr=nsig.true.fdr>prctile.nsig.fdr[1]
  
  
  # test other stat: orignal r, r sqr, r abs
  cor.true.orig=cor.true
  cor.true.sqr=cor.true**2
  cor.true.abs=abs(cor.true)
  
  stat.true.orig=mean(cor.true.orig)
  stat.true.sqr=mean(cor.true.sqr)
  stat.true.abs=mean(cor.true.abs)
  
  
  cor.null.orig=cor.null
  cor.null.sqr=cor.null**2
  cor.null.abs=abs(cor.null)
  
  stat.null.orig=apply(cor.null.orig,2,mean)
  stat.null.sqr=apply(cor.null.sqr,2,mean)
  stat.null.abs=apply(cor.null.abs,2,mean)
  
  prctile.orig=quantile(stat.null.orig, c(0.025, 0.975))#lower 2.5% and upper 97.5%
  prctile.sqr=quantile(stat.null.sqr, c(0.95))
  prctile.abs=quantile(stat.null.abs, c(0.95))
  
  if_sig.orig=!(prctile.orig[1]<stat.true.orig && stat.true.orig<prctile.orig[2]) # not inside 95%
  if_sig.sqr=stat.true.sqr>=prctile.sqr[1] # greater than 95%
  if_sig.abs=stat.true.abs>=prctile.abs[1] # greater than 95%
  
  
  mean.under.null.orig=mean(stat.null.orig)
  #mean.under.null.sqr=mean(stat.null.sqr)
  #mean.under.null.abs=mean(stat.null.abs)
  
  stat.null.oirg.aug=c(stat.true.orig,stat.null.orig)
  stat.null.sqr.aug=c(stat.true.sqr,stat.null.sqr)
  stat.null.abs.aug=c(stat.true.abs,stat.null.abs)
  
  stat.true.orig.adj=stat.true.orig-mean.under.null.orig
  #stat.true.sqr.adj=stat.true.sqr-mean.under.null.sqr
  #stat.true.abs.adj=stat.true.abs-mean.under.null.abs
  
  
  stat.null.orig.aug.adj=stat.null.oirg.aug-mean.under.null.orig
  #stat.null.sqr.aug.adj=stat.null.sqr.aug-mean.under.null.sqr
  #stat.null.abs.aug.adj=stat.null.abs.aug-mean.under.null.abs
  
  p.orig=sum(abs(stat.null.orig.aug.adj) >=abs(stat.true.orig.adj))/length(stat.null.orig.aug.adj)
  p.sqr=sum(stat.null.sqr.aug>=stat.true.sqr)/length(stat.null.sqr.aug)
  p.abs=sum(stat.null.abs.aug>=stat.true.abs)/length(stat.null.abs.aug)
  
  
  res_df=data.frame(mean_r=stat.true.orig,
                    mean_r.null025=as.numeric(prctile.orig[1]),
                    mean_r.null975=as.numeric(prctile.orig[2]),
                    mean_r.if_sig=as.logical(if_sig.orig),
                    mean_r.adj=stat.true.orig.adj,
                    mean_r.p=p.orig,
                    
                    mean_rsqr=stat.true.sqr,
                    mean_rsqr.null095=as.numeric(prctile.sqr[1]),
                    mean_rsqr.if_sig=as.logical(if_sig.sqr),
                    mean_rsqr.p=p.sqr,
                    
                    mean_rabs=stat.true.abs,
                    mean_rabs.null095=as.numeric(prctile.abs[1]),
                    mean_rabs.if_sig=as.logical(if_sig.abs),
                    mean_rabs.p=p.abs,
                    
                    nsig=as.numeric(nsig.true),
                    nsig.null095=as.numeric(prctile.nsig[1]),
                    nsig.if_sig=as.logical(if_sig.nsig),
                    nsig.fdr=as.numeric(nsig.true.fdr),
                    nsig.fdr.null095=as.numeric(prctile.nsig.fdr[1]),
                    nsig.if_sig.fdr=as.logical(if_sig.nsig.fdr))
  return (res_df)
}

caculate_FPR <- function(glist_expr_df,
                         gall_expr_df,
                         random_brain4FPR,
                         shuffle.id=NULL){
  ntest=dim(random_brain4FPR)[2]
  all_if_sig.orig=list()
  all_if_sig.sqr=list()
  all_if_sig.abs=list()
  all_if_sig.nsig=list()
  all_if_sig.nsig.fdr=list()
  for (test_i in c(1:ntest)){
    random_brain_df=data.frame(x=random_brain4FPR[,test_i],
                               row.names = rownames(gall_expr_df))
    
    res_df=sigtest(random_brain_df,
                    glist_expr_df,
                    gall_expr_df,
                    shuffle.id)
    all_if_sig.orig[[paste0('ntest_',test_i)]] = res_df$mean_r.if_sig
    all_if_sig.sqr[[paste0('ntest_',test_i)]] = res_df$mean_rsqr.if_sig
    all_if_sig.abs[[paste0('ntest_',test_i)]] = res_df$mean_rabs.if_sig
    all_if_sig.nsig[[paste0('ntest_',test_i)]] = res_df$nsig.if_sig
    all_if_sig.nsig.fdr[[paste0('ntest_',test_i)]] = res_df$nsig.if_sig.fdr
  }
  FPR.orig=Reduce(sum,all_if_sig.orig)/ntest
  FPR.sqr=Reduce(sum,all_if_sig.sqr)/ntest
  FPR.abs=Reduce(sum,all_if_sig.abs)/ntest
  FPR.nsig=Reduce(sum,all_if_sig.nsig)/ntest
  FPR.nsig.fdr=Reduce(sum,all_if_sig.nsig.fdr)/ntest
  return(c(FPR.orig,FPR.sqr,FPR.abs,FPR.nsig,FPR.nsig.fdr))
}




simulating_glist<- function(glist_n=5,
                            sim_rep=100,
                            gall_expr_df,
                            SIMseed=2022,
                            random_brain4FPR,
                            shuffle.id=NULL){
  # SIMseed should change for each iteration if run a for loop for glist_n with sim_rep=1
  sim_res=list()
  set.seed(SIMseed)
  seeds=sample(100000:999999,size=sim_rep,replace = F)
  for (sim_i in c(1:sim_rep)){
    cat(sprintf('Sim_%s/%s\n',sim_i,sim_rep))
    seed2use=seeds[sim_i]
    set.seed(seed2use)
    glist_idx=sample(dim(gall_expr_df)[2],glist_n,replace=F)
    glist=colnames(gall_expr_df)[glist_idx]
    glist_expr_df=get_glist_expression(gall_expr_df,glist)
    sim_res[['seed']][sim_i]=seed2use
    sim_res[['glist_n']][sim_i]=glist_n

    #sim_res[['glist_moran']][sim_i]=colMeans(get_df_moran_idx(glist_expr_df))

    FPR_vect=caculate_FPR(glist_expr_df,
                          gall_expr_df, 
                          random_brain4FPR, 
                          shuffle.id=shuffle.id)
    
    sim_res[['FPR.orig']][sim_i]=FPR_vect[1]
    sim_res[['FPR.sqr']][sim_i]=FPR_vect[2]
    sim_res[['FPR.abs']][sim_i]=FPR_vect[3]
    sim_res[['FPR.nsig']][sim_i]=FPR_vect[4]
    sim_res[['FPR.nsig.fdr']][sim_i]=FPR_vect[5]
    
    coexp_df=caculate_coexp(glist_expr_df)
    sim_res[['coexp']][sim_i]=coexp_df$mean
    sim_res[['coexp_pos_adj']][sim_i]=coexp_df$adj_posmean
    sim_res[['coexp_neg_adj']][sim_i]=coexp_df$adj_negmean
    
  }
  sim_res_df=data.frame(sim_res)
  return(sim_res_df)
}    




## functions for synGO
# test real GO term

split_gene <- function(gene_string,pattern=';'){
  gene_list=unlist(strsplit(gene_string, pattern))
  return(gene_list)}

combine_gene <- function(gene_list,pattern=';'){
  gene_string=paste0(gene_list,collapse = pattern)
  return(gene_string)
}



find_gene_n <- function(gene_string,pattern=';'){
  #gene_string=as.character(gene_string)
  gene_list=unlist(strsplit(gene_string, pattern))
  gene_n=length(gene_list)
  return(gene_n)
}



find_gene_in_AHBA <- function(gene_string, 
                              AHBA_glist,
                              pattern=';'){
  gene_list=unlist(strsplit(gene_string, pattern))
  int_gene=intersect(gene_list, AHBA_glist)
  int_gene_string=combine_gene(int_gene)
  return(int_gene_string)
}


# # Function to generate a permutation map from a set of cortical regions of interest to itself, 
# # while (approximately) preserving contiguity and hemispheric symmetry.
# # The function is based on a rotation of the FreeSurfer projection of coordinates
# # of a set of regions of interest on the sphere.
# #
# # Inputs:
# # coord.l       coordinates of left hemisphere regions on the sphere        array of size [n(LH regions) x 3]
# # coord.r       coordinates of right hemisphere regions on the sphere       array of size [n(RH regions) x 3]
# # nrot          number of rotations (default = 10000)                       scalar
# #
# # Output:
# # perm.id      array of permutations, from set of regions to itself        array of size [n(total regions) x nrot]
# #
# # required library: matrixStats (for rowMins function)
# #
# # Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
# #       Updated on 16/10/2019 with permutation scheme that uniformly samples the space of permutations on the sphere
# #       See github repo (@frantisekvasa) and references within for details
# 
# library(matrixStats)
# rotate.parcellation.lh = function(coord.l,nrot=10000) {
#   # check that coordinate dimensions are correct
#   if (!all(dim(coord.l)[2]==3)) {
#     if (all(dim(coord.l)[1]==3)) {
#       print('transposing coordinates to be of dimension nROI x 3')
#       coord.l = t(coord.l)
#     }
#   }
#   
#   nroi.l = dim(coord.l)[1]   # n(regions) in the left hemisphere
#   nroi = nroi.l      # total n(regions)
#   
#   perm.id = array(0,dim=c(nroi,nrot)); # initialise output array
#   r = 0; c = 0; # count successful (r) and unsuccessful (c) iterations
#   
#   # UPDATED 16/10/2019 - set up updated permutation scheme 
#   I1 = diag(3); I1[1,1] = -1;
#   # main loop -  use of "while" is to ensure any rotation that maps to itself is excluded (this is rare, but can happen)
#   while (r < nrot) {
#     
#     # UPDATED 16/10/2019
#     A = matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
#     qrdec = qr(A)       # QR decomposition
#     TL = qr.Q(qrdec)    # Q matrix
#     temp = qr.R(qrdec)  # R matrix
#     TL = TL%*%diag(sign(diag(temp)))
#     if (det(TL)<0) {
#       TL[,1] = -TL[,1]
#     }
#     # reflect across the Y-Z plane for right hemisphere
#     
#     coord.l.rot = coord.l %*% TL; # transformed (rotated) left coordinates
#     dist.l = array(0,dim=c(nroi.l,nroi.l));
#     for (i in 1:nroi.l) { # left
#       for (j in 1:nroi.l) {
#         dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
#       }
#     }
#     
#     
#     # LEFT
#     # calculate distances, proceed in order of "most distant minimum"
#     # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
#     # as this region is the hardest to match and would only become harder as other regions are assigned
#     temp.dist.l = dist.l
#     rot.l = c(); ref.l = c();
#     #tba.r = tba.c = 1:nroi.l # rows and columns that are yet "to be assigned"
#     for (i in 1:nroi.l) {
#       # max(min) (described above)
#       ref.ix = which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   # "furthest" row
#       rot.ix = which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) ) # closest region
#       
#       # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
#       # ref.ix = which(nanmean(temp.dist.l,2)==nanmax(nanmean(temp.dist.l,2)))    # "furthest" row
#       # rot.ix = which(temp.dist.l(ref.ix,:)==nanmin(temp.dist.l(ref.ix,:)))      # closest region    
#       ref.l = c(ref.l,ref.ix) # store reference and rotated indices
#       rot.l = c(rot.l,rot.ix)
#       temp.dist.l[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
#       temp.dist.l[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
#       #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
#       #temp.dist.l[ref.ix,] = NA
#     }
#     
#     
#     # mapping is x->y
#     # collate vectors from both hemispheres + sort mapping according to "reference" vector
#     b = sort(ref.l,index.return=T); 
#     ref.l.sort = ref.l[b$ix]; rot.l.sort = rot.l[b$ix];
#     
#     # verify that permutation worked (output should be vector with values 1:nroi = 1:(nroi_l+nroi_r))
#     if (!all(sort(rot.l.sort,decreasing=F)==c(1:nroi))) {
#       #save.image('~/Desktop/perm_error.RData')
#       browser("permutation error")
#     }
#     
#     # verify that permutation does not map to itself
#     if (!all(rot.l.sort==c(1:nroi))) {
#       r = r+1
#       perm.id[,r] = rot.l.sort # if it doesn't, store it
#     } else {
#       c = c+1
#       print(paste('map to itself n. ',toString(c),sep=''))
#     }
#     
#     # track progress
#     if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
#   }
#   return(perm.id)
# }



# get_34perm.id <- function(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data',
#                           nreps=10000,
#                           if_spin=T){
#   if (if_spin){
#     lh_centroid_data=read.csv(sprintf('%s/lh_centroid.csv',data_path))
#     centroid.l=lh_centroid_data %>%  select(starts_with('centroid'))
#     perm.id=rotate.parcellation.lh(as.matrix(centroid.l),
#                                    nrot=nreps+100) # take 100 more in case there are duplicated cols
#   }else {
#     id.mat=matrix(rep(c(1:34),times=nreps+100),nrow = 34,ncol = nreps+100)
#     perm.id=apply(id.mat, MARGIN = 2, function(x) sample(x, replace = FALSE, size = length(x)))}
#   
#   # B = matrix(c(1, 4, 0, 2, 56, 7, 1, 4, 0, 33, 2, 5), nrow = 3)
#   # colnames(B) <- c("a", "b", "c", "d")
#   # 
#   # B
#   # ##      a  b c  d
#   # ## [1,] 1  2 1 33
#   # ## [2,] 4 56 4  2
#   # ## [3,] 0  7 0  5
#   # 
#   # B[, !duplicated(t(B))]
#   # ##      a  b  d
#   # ## [1,] 1  2 33
#   # ## [2,] 4 56  2
#   # ## [3,] 0  7  5
#   print(paste0('find duplciated columns: ', sum(duplicated(t(perm.id)))))
#   perm.id.nodup=perm.id[, !duplicated(t(perm.id))]
#   perm.id.final=perm.id.nodup[,1:nreps]
#   return(perm.id.final)
# }

# create perm id and save on disk
# perm.id.spin10k=get_34perm.id(nreps = 10000, if_spin = T)
# perm.id.nospin10k=get_34perm.id(nreps = 10000, if_spin =F)
# saveRDS(perm.id.spin10k,paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_spin10k.rds',sep='/'))
# saveRDS(perm.id.nospin10k,paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_nospin10k.rds',sep='/'))

# Read
# perm.id.spin10k=readRDS(paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_spin10k.rds',sep='/')) 
# perm.id.nospin10k=readRDS(paste('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data','34perm_id_nospin10k.rds',sep='/')) 




rbind_sim_results<-function(sim_path){
  files.names<-list.files(sim_path,pattern="*.csv")
  print(files.names)
  # read csv files and return them as items in a list()
  theList <- lapply(files.names, function(x){read.csv(sprintf('%s/%s',sim_path,x))})
  theResult <- do.call(rbind,theList)
  write.csv(theResult, sprintf('%s/all_sim_results.csv',sim_path), row.names = F)
}


rbind_sim_results1<-function(sim_path,filepattern){
  pattern2search=sprintf('^sim_%s_.*\\.csv$',filepattern)
  files.names<-list.files(sim_path,pattern=pattern2search)
  print(files.names)
  files.n=length(files.names)
  # read csv files and return them as items in a list()
  theList <- lapply(files.names, function(x){read.csv(sprintf('%s/%s',sim_path,x))})
  theResult <- do.call(rbind,theList)
  write.csv(theResult, sprintf('%s/all_%s_sims_%s.csv',sim_path,files.n,filepattern), row.names = F)
}

