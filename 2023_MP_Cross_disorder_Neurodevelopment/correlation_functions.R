# correlation functions



get_68perm.id <- function(nrot=10000,seed=666){
  lh_centroid_data=read.csv('F:/Google Drive/post-doc/p_factor/summary_stat/allenbrain/spin_test/lh_centroid.csv')
  rh_centroid_data=read.csv('F:/Google Drive/post-doc/p_factor/summary_stat/allenbrain/spin_test/rh_centroid.csv')
  lr_centroid=rbind(lh_centroid_data,rh_centroid_data)
  colnames(lr_centroid)=c('label','centroid1','centroid2','centroid3')
  centroid.l=lh_centroid_data %>%  select(starts_with('centroid'))
  centroid.r=rh_centroid_data %>%  select(starts_with('centroid'))  
  perm.id=rotate.parcellation(as.matrix(centroid.l),as.matrix(centroid.r),nrot=nrot)
  
  return(perm.id)
}


caculate_spin_p <- function(df2cor,perm.id){
  p.spin <-matrix(nrow = ncol(df2cor), ncol = ncol(df2cor))
  for (i in (1:ncol(df2cor))) {
    for (j in (i:ncol(df2cor))){
      if (i==j){p.spin[i,j]=NA}
      else{
        p.spin[i,j] <- perm.sphere.p(df2cor[,i],df2cor[,j],perm.id,corr.type='pearson')}
    }
  }
  #       [,1] [,2] [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
  # [1,]    0    0    0 0.005 0.005 0.120 0.055 0.000 0.010 0.000
  # [2,]   NA    0    0 0.000 0.005 0.485 0.000 0.000 0.010 0.000
  # [3,]   NA   NA    0 0.005 0.000 0.180 0.010 0.000 0.000 0.000
  # [4,]   NA   NA   NA 0.000 0.005 0.000 0.555 0.150 0.455 0.255
  # [5,]   NA   NA   NA    NA 0.000 0.320 0.025 0.000 0.000 0.005
  # [6,]   NA   NA   NA    NA    NA 0.000 0.090 0.015 0.000 0.000
  # [7,]   NA   NA   NA    NA    NA    NA 0.000 0.190 0.005 0.030
  # [8,]   NA   NA   NA    NA    NA    NA    NA 0.000 0.000 0.000
  # [9,]   NA   NA   NA    NA    NA    NA    NA    NA 0.000 0.000
  # [10,]   NA   NA   NA    NA    NA    NA    NA    NA    NA 0.000
  
  
  # copy upper to lower
  p.spin[lower.tri(p.spin)]=t(p.spin)[lower.tri(p.spin)]
  #       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
  # [1,] 0.000 0.000 0.000 0.005 0.005 0.120 0.055 0.000 0.010 0.000
  # [2,] 0.000 0.000 0.000 0.000 0.005 0.485 0.000 0.000 0.010 0.000
  # [3,] 0.000 0.000 0.000 0.005 0.000 0.180 0.010 0.000 0.000 0.000
  # [4,] 0.005 0.000 0.005 0.000 0.005 0.000 0.555 0.150 0.455 0.255
  # [5,] 0.005 0.005 0.000 0.005 0.000 0.320 0.025 0.000 0.000 0.005
  # [6,] 0.120 0.485 0.180 0.000 0.320 0.000 0.090 0.015 0.000 0.000
  # [7,] 0.055 0.000 0.010 0.555 0.025 0.090 0.000 0.190 0.005 0.030
  # [8,] 0.000 0.000 0.000 0.150 0.000 0.015 0.190 0.000 0.000 0.000
  # [9,] 0.010 0.010 0.000 0.455 0.000 0.000 0.005 0.000 0.000 0.000
  # [10,] 0.000 0.000 0.000 0.255 0.005 0.000 0.030 0.000 0.000 0.000
  
  nam=dimnames(df2cor)[[2]]
  dimnames(p.spin) =list(nam,nam)
  return(p.spin)
}


 # plot correlation matrix and show signficant correlation with spin_p<0.05
cor_plot <- function(cor.list,...){
  # type: upper or lower 
  col2use=colorRampPalette(c("steelblue1","white","firebrick1"))
  corrplot(cor.list$r, method="color", col=col2use(200),  
           addCoef.col = "black", 
           tl.col="black", tl.srt=45, 
           p.mat = cor.list$p.spin, sig.level = 0.05, insig = "blank", 
           diag=FALSE ,
           number.cex=0.8,cl.ratio = 0.3, addgrid.col='grey',...
  )
}


# correlation results to save
cor_t2save <- function(cor.list){
  r2save=cor.list$r
  p2save=cor.list$P
  p.spin2save=cor.list$p.spin
  
  r2save[lower.tri(r2save)] <- NA
  p2save[lower.tri(p2save)] <- NA
  p.spin2save[lower.tri(p.spin2save)] <- NA
  t2save=rbind(r2save,p2save,p.spin2save)
  return(t2save)
}


# shuffle loading according to perm.id
build_null_brain_df <- function(perm.id,brain_df){
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


# p is calculated as the quantile of the true r in the null distribution
# note: a better way is to +1 so that p will not be 0.

caculate_p_from_null <- function(pos_null_dist,neg_null_dist, r_val){
  if (r_val>0){
    p_val=sum(pos_null_dist>r_val)/length(pos_null_dist)
  }else{
    p_val=sum(neg_null_dist<r_val)/length(neg_null_dist)
  }
  return(p_val)
}

