# external functions

# syndromics functions
stand_loadings<-function(pca,pca_data){
  if (class(pca)[1]=='prcomp'){
    if (is.numeric(pca$scale)){
      loadings<-as.data.frame((pca$rotation%*%diag(pca$sdev)))
    }else{
      loadings<-as.data.frame((pca$rotation%*%diag(pca$sdev))/apply(pca_data,2,sd))
    }
    colnames(loadings)<-paste('PC', 1:dim(pca$x)[2],sep = '')
  }else if(class(pca)[1]%in%"princals"){
    loadings<-as.data.frame(pca$loadings)
    colnames(loadings)<-paste('PC', 1:pca$ndim,sep = '')
  }else{
    stop('pca must be a prcomp object')
  }
  return(loadings)
}

VAF_plot<-function(pca, ndim=1:5, resample_ci=NULL, style="line", colors=c("steelblue","orange")){

  if (class(pca)[1]=='prcomp'){
    original_VAF<-pca$sdev^2/sum(pca$sdev^2)
  }else if (class(pca)[1]%in%c('princals')){
    original_VAF<-pca$evals/sum(pca$evals)
  }

  plot_df<-data.frame("value"=original_VAF[ndim],
                      "component"=factor(paste0('PC',ndim), levels = paste0('PC',ndim)),
                      "name"="Original")

  if (style=="line"){

    vaf_plot<-ggplot(plot_df, aes(.data$component,.data$value*100))+
      geom_point(aes(color="Original"), show.legend = F)+
      geom_line(aes(group=.data$name, color="Original"), show.legend = F)+
      theme_minimal()+
      xlab(NULL)+ylab('Variance accounted for (VAF %)')+
      theme(legend.position = c(0.7,0.8), legend.background = element_rect(color='grey'))+
      scale_colour_manual(values=colors[1])

    if (!is.null(resample_ci)){
      resample_ci<-resample_ci%>%
        mutate(component=factor(rownames(resample_ci), levels = rownames(resample_ci)),
               name="Permuted")
      resample_ci<-resample_ci[ndim,]
      suppressMessages(vaf_plot<-vaf_plot+
        geom_point(data=resample_ci, aes(x=.data$component, y=.data$mean*100, color="Permuted"),
                   inherit.aes = F)+
        geom_line(data=resample_ci, aes(x=.data$component, y=.data$mean*100,group=.data$name,color="Permuted"),
                  inherit.aes = F, show.legend = TRUE)+
        geom_errorbar(data=resample_ci,
                      aes(x=.data$component, y=.data$mean*100, ymin=.data$ci_low*100,ymax=.data$ci_high*100,color="Permuted"),
                      inherit.aes = F,width=0.5)+
        scale_color_manual(values=c(colors[1],colors[2]))+
        theme(legend.title = element_blank()))

    }
  }
  if(style=="reduced"){
    plot_df$VAF_value<-
      factor(paste0(round(plot_df$value*100,1),"%"),
             levels=rev(paste0(round(plot_df$value*100,1),"%")))

    vaf_plot<-ggplot(plot_df, aes(.data$value, .data$VAF_value))+
      geom_segment(aes(xend=0, yend=.data$VAF_value), size=1.5, color=colors[1])+
      xlab(NULL)+ylab('VAF (%)')+
      scale_x_continuous(breaks = NULL)+
      theme_minimal()+
      theme(panel.grid = element_blank())

    if(!is.null(resample_ci)){
      resample_ci$component<-rownames(resample_ci)
      resample_ci$name<-"Permuted"
      resample_ci$VAF_value<-
        factor(paste0(round(plot_df$value*100,1),"%"),
               levels=rev(paste0(round(plot_df$value*100,1),"%")))
      vaf_plot<-vaf_plot+
        geom_errorbar(data = resample_ci, aes(xmin=.data$ci_low, xmax=.data$ci_high,
                                              x=.data$mean, y=.data$VAF_value),
                      width=0.4,inherit.aes = F, color=colors[2])+
        geom_point(data = resample_ci, aes(x=.data$mean, y=.data$VAF_value),
                   inherit.aes = F, color=colors[2])
    }
  }
  return(vaf_plot)
}




barmap_loading_rename<-function(pca.obj, ndim=1:3, cutoff=0,resample_ci=NULL,conf=0.95, plot_list_center=F,
                         plot_legend= TRUE,text_values=F, star_values=F,
                         text_size=4, plot_cutoff= TRUE, vars=NULL, colors=c("steelblue1","white","firebrick1"),
                         gradient_color= TRUE){
  

  load_df=pca.obj$renamed_loading
  load_df$Variables=rownames(load_df)
  old_scipen<-getOption('scipen')
  on.exit(options(scipen=old_scipen))
  
  options(scipen=999)
  
  vars=rev(rownames(load_df))# reverse order
   
  if (max(ndim)>ncol(load_df[,-1])){
    stop('higher value for ndim can not be higher than the number of PCs')
  }
  
  if (!is.null(vars)){
    if(!is.character(vars)){
      stop('vars must be a character vector')
    }
   
    order_var<-match(vars, load_df$Variables)
    if (sum(is.na(order_var))>0){
      stop(paste('variable/s not found: ',
                 vars[which(is.na(order_var))], sep = ''))
    }
    
    load_df<-load_df%>%filter(.data$Variables%in%vars)%>%
      mutate(Variables=factor(.data$Variables,vars))
  }
  
  if (length(cutoff)!=1 & length(cutoff)!=length(ndim)){
    stop('cutoff length must be one or equal to ndim length')
  }
  
  cutoff_df<-load_df%>%pivot_longer(-(.data$Variables),
                                    names_to = "component", values_to = "loading")%>%
    filter(.data$component%in%paste('PC',ndim, sep = ''))%>%
    group_by(.data$component)%>%
    summarise(count=n())%>%mutate(cutoff=cutoff)%>%
    mutate(component=factor(.data$component, levels=names(load_df)))
  
  load_df<-load_df%>%pivot_longer(-(.data$Variables),
                                  names_to = "component", values_to = "loading")%>%
    filter(.data$component%in%paste('PC',ndim, sep = ''))%>%
    left_join(cutoff_df, by='component')%>%
    group_by(.data$component)%>%
    mutate(cutoff=cutoff,
           star=ifelse(abs(.data$loading)>=cutoff,'*',''),
           loading_txt=as.character(round(.data$loading,2)),
           weight=(round(.data$loading,2)))%>%
    ungroup()%>%
    mutate(component=factor(.data$component, levels=names(load_df)))
  
  b_plot<-load_df%>%
    ggplot(aes(.data$Variables, .data$loading))+
    geom_col(show.legend = plot_legend)+
    theme_minimal()+
    coord_flip()+
    facet_grid(~.data$component)
  
  #PCA list
  if (!is.null(resample_ci)){
    
    if (max(ndim)>length(unique(resample_ci$component))){
      stop(paste0('higher value for ndim can not be higher than the
                  number of PCs in resample_ci which is ',length(unique(resample_ci$component))))
    }
    
    resample_ci<-as.data.frame(resample_ci%>%
                                 filter(.data$component%in%unique(load_df$component)))
    
    b_plot<-b_plot+
      geom_errorbar(data=resample_ci, aes(x=.data$Variables,ymin=.data$ci_low, ymax=.data$ci_high),width=0.5, inherit.aes = F)
    
    if (plot_list_center){
      b_plot<-b_plot+
        geom_point(data=resample_ci, aes(y=.data$mean, x=.data$Variables),inherit.aes = F)
    }
    }
  
  if (text_values){
    b_plot<-b_plot+geom_text(aes(label=.data$loading_txt,y=.data$loading*1.1),size=2)
  }else if (star_values){
    b_plot<-b_plot+geom_text(aes(label=.data$star, y=.data$loading*1.1), color='black')
  }
  
  if (gradient_color){
    b_plot<-b_plot+aes(.data$Variables, .data$loading, fill=.data$loading)+
      scale_fill_gradient2(low=colors[1], high=colors[3],mid = colors[2],
                           limits=c(-1,1), na.value =  "transparent")
  }else{
    b_plot<-b_plot+aes(.data$Variables, .data$loading, fill=(.data$loading)>0)+
      scale_fill_manual(values=c(colors[1],colors[3]), na.value =  "transparent")+
      theme(legend.position = "none")
  }
  
  if(plot_cutoff){
    b_plot<-b_plot+
      geom_hline(data=cutoff_df,aes(yintercept = cutoff), color=colors[3], alpha=0.6)+
      geom_hline(data=cutoff_df,aes(yintercept = -cutoff), color=colors[1], alpha=0.6)
  }
  
  b_plot<-b_plot+ylim(-1.1,1.1)+
    labs(fill='s.loading')+
    geom_hline(yintercept = 0, color='black')+
    ylab(NULL)+xlab(NULL)
  
  return(b_plot)
}











# Function to generate a permutation map from a set of cortical regions of interest to itself, 
# while (approximately) preserving contiguity and hemispheric symmetry.
# The function is based on a rotation of the FreeSurfer projection of coordinates
# of a set of regions of interest on the sphere.
#
# Inputs:
# coord.l       coordinates of left hemisphere regions on the sphere        array of size [n(LH regions) x 3]
# coord.r       coordinates of right hemisphere regions on the sphere       array of size [n(RH regions) x 3]
# nrot          number of rotations (default = 10000)                       scalar
#
# Output:
# perm.id      array of permutations, from set of regions to itself        array of size [n(total regions) x nrot]
#
# required library: matrixStats (for rowMins function)
#
# Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
#       Updated on 16/10/2019 with permutation scheme that uniformly samples the space of permutations on the sphere
#       See github repo (@frantisekvasa) and references within for details

rotate.parcellation = function(coord.l,coord.r,nrot=10000) {
  
  # check that coordinate dimensions are correct
  if (!all(dim(coord.l)[2]==3,dim(coord.r)[2]==3)) {
    if (all(dim(coord.l)[1]==3,dim(coord.r)[1]==3)) {
      print('transposing coordinates to be of dimension nROI x 3')
      coord.l = t(coord.l)
      coord.r = t(coord.r)
    }
  }
  
  nroi.l = dim(coord.l)[1]   # n(regions) in the left hemisphere
  nroi.r = dim(coord.r)[1]   # n(regions) in the right hemisphere
  nroi = nroi.l+nroi.r       # total n(regions)
  
  perm.id = array(0,dim=c(nroi,nrot)); # initialise output array
  r = 0; c = 0; # count successful (r) and unsuccessful (c) iterations
  
  # UPDATED 16/10/2019 - set up updated permutation scheme 
  I1 = diag(3); I1[1,1] = -1;
  # main loop -  use of "while" is to ensure any rotation that maps to itself is excluded (this is rare, but can happen)
  while (r < nrot) {
    
    # UPDATED 16/10/2019
    A = matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec = qr(A)       # QR decomposition
    TL = qr.Q(qrdec)    # Q matrix
    temp = qr.R(qrdec)  # R matrix
    TL = TL%*%diag(sign(diag(temp)))
    if (det(TL)<0) {
      TL[,1] = -TL[,1]
    }
    # reflect across the Y-Z plane for right hemisphere
    TR = I1 %*% TL %*% I1;
    coord.l.rot = coord.l %*% TL; # transformed (rotated) left coordinates
    coord.r.rot = coord.r %*% TR; # transformed (rotated) right coordinates
    
    # OLD PERMUTATION SCHEME - COMMENTED OUT 16/10/2019
    # # choose three angles at random - x, y, z
    # # random angles
    # rm(.Random.seed, envir=globalenv()) # reset random seed
    # ax = 2*pi*runif(1)
    # ay = 2*pi*runif(1)
    # az = 2*pi*runif(1)
    # 
    # ### rotation matrices
    # # left hemisphere
    # rx.l = cbind(c(1,0,0),c(0,cos(ax),sin(ax)),c(0,-sin(ax),cos(ax)))
    # ry.l = cbind(c(cos(ay),0,-sin(ay)),c(0,1,0),c(sin(ay),0,cos(ay)))
    # rz.l = cbind(c(cos(az),sin(az),0),c(-sin(az),cos(az),0),c(0,0,1))
    # # right hemisphere - same magnitude of rotation, but signs for y and z axes are flipped to retain symmetry
    # rx.r = cbind(c(1,0,0),c(0,cos(ax),sin(ax)),c(0,-sin(ax),cos(ax)))
    # ry.r = cbind(c(cos(-ay),0,-sin(-ay)),c(0,1,0),c(sin(-ay),0,cos(-ay)))
    # rz.r = cbind(c(cos(-az),sin(-az),0),c(-sin(-az),cos(-az),0),c(0,0,1))
    # 
    # # perform rotation (mutiply coordinates by rotation matrices, for n-by-3 matrix)
    # # left hemisphere
    # coord.l.rot.x = coord.l %*% rx.l
    # coord.l.rot.xy = coord.l.rot.x %*% ry.l
    # coord.l.rot.xyz = coord.l.rot.xy %*% rz.l
    # # right hemisphere
    # coord.r.rot.x = coord.r %*% rx.r
    # coord.r.rot.xy = coord.r.rot.x %*% ry.r
    # coord.r.rot.xyz = coord.r.rot.xy %*% rz.r
    
    # after rotation, find "best" match between rotated and unrotated coordinates
    # first, calculate distance between initial coordinates and rotated ones
    dist.l = array(0,dim=c(nroi.l,nroi.l));
    dist.r = array(0,dim=c(nroi.r,nroi.r));
    # OLD PERMUTATION SCHEME - COMMENTED OUT 16/10/2019
    # for (i in 1:nroi.l) { # left
    #   for (j in 1:nroi.l) {
    #     dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot.xyz[j,])^2 ) )
    #   }
    # }
    # for (i in 1:nroi.r) { # right
    #   for (j in 1:nroi.r) {
    #     dist.r[i,j] = sqrt( sum( (coord.r[i,]-coord.r.rot.xyz[j,])^2 ) )
    #   }
    # }
    # UPDATED 5/9/2019 - change of rotated variable name to "coord.l/r.rot" (from coord.l/r.rot.xyz)
    for (i in 1:nroi.l) { # left
      for (j in 1:nroi.l) {
        dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
    for (i in 1:nroi.r) { # right
      for (j in 1:nroi.r) {
        dist.r[i,j] = sqrt( sum( (coord.r[i,]-coord.r.rot[j,])^2 ) )
      }
    }
    
    # LEFT
    # calculate distances, proceed in order of "most distant minimum"
    # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
    # as this region is the hardest to match and would only become harder as other regions are assigned
    temp.dist.l = dist.l
    rot.l = c(); ref.l = c();
    #tba.r = tba.c = 1:nroi.l # rows and columns that are yet "to be assigned"
    for (i in 1:nroi.l) {
      # max(min) (described above)
      ref.ix = which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   # "furthest" row
      rot.ix = which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) ) # closest region
      
      # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
      # ref.ix = which(nanmean(temp.dist.l,2)==nanmax(nanmean(temp.dist.l,2)))    # "furthest" row
      # rot.ix = which(temp.dist.l(ref.ix,:)==nanmin(temp.dist.l(ref.ix,:)))      # closest region    
      ref.l = c(ref.l,ref.ix) # store reference and rotated indices
      rot.l = c(rot.l,rot.ix)
      temp.dist.l[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
      temp.dist.l[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
      #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
      #temp.dist.l[ref.ix,] = NA
    }
    
    # RIGHT
    # calculate distances, proceed in order of "most distant minimum"
    # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
    # as this region is the hardest to match and would only become harder as other regions are assigned
    temp.dist.r = dist.r;
    rot.r = c(); ref.r = c();
    for (i in 1:nroi.r) {
      # max(min) (described above)
      ref.ix = which( rowMins(temp.dist.r,na.rm=T) == max(rowMins(temp.dist.r,na.rm=T),na.rm=T) )   # "furthest" row
      rot.ix = which( temp.dist.r[ref.ix,] == min(temp.dist.r[ref.ix,],na.rm=T) )             # closest region
      
      # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
      # ref.ix = which(nanmean(temp.dist.r,2)==nanmax(nanmean(temp.dist.r,2)))    # "furthest" row
      # rot.ix = which(temp.dist.r(ref.ix,:)==nanmin(temp.dist.r(ref.ix,:)))      # closest region
      ref.r = c(ref.r,ref.ix) # store reference and rotated indices
      rot.r = c(rot.r,rot.ix)
      temp.dist.r[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
      temp.dist.r[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
      #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
      #temp.dist.l[ref.ix,] = NA
    }
    
    # mapping is x->y
    # collate vectors from both hemispheres + sort mapping according to "reference" vector
    ref.lr = c(ref.l,nroi.l+ref.r); rot.lr = c(rot.l,nroi.l+rot.r);
    b = sort(ref.lr,index.return=T); 
    ref.lr.sort = ref.lr[b$ix]; rot.lr.sort = rot.lr[b$ix];
    
    # verify that permutation worked (output should be vector with values 1:nroi = 1:(nroi_l+nroi_r))
    if (!all(sort(rot.lr.sort,decreasing=F)==c(1:nroi))) {
      #save.image('~/Desktop/perm_error.RData')
      browser("permutation error")
    }
    
    # verify that permutation does not map to itself
    if (!all(rot.lr.sort==c(1:nroi))) {
      r = r+1
      perm.id[,r] = rot.lr.sort # if it doesn't, store it
    } else {
      c = c+1
      print(paste('map to itself n. ',toString(c),sep=''))
    }
    
    # track progress
    if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
    
  }
  
  return(perm.id)
  
}




# Function to generate a p-value for the spatial correlation between two parcellated cortical surface maps, 
# using a set of spherical permutations of regions of interest (which can be generated using the function "rotate_parcellation").
# The function performs the permutation in both directions; i.e.: by permute both measures, 
# before correlating each permuted measure to the unpermuted version of the other measure
#
# Inputs:
# x                 one of two maps to be correlated                                                                    vector
# y                 second of two maps to be correlated                                                                 vector
# perm_id           array of permutations, from set of regions to itself (as generated by "rotate_parcellation")        array of size [n(total regions) x nrot]
# corr_type         type of correlation                                                                                 "spearman" (default) or "pearson"
#
# Output:
# p_perm            permutation p-value
#
# Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018

perm.sphere.p = function(x,y,perm.id,corr.type='spearman') {
  
  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations
  
  rho.emp = cor(x,y,method=corr.type)  # empirical correlation
  
  # permutation of measures
  x.perm = y.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] = x[perm.id[i,r]]
      y.perm[i,r] = y[perm.id[i,r]]
    }
  }
  
  # correlation to unpermuted measures
  rho.null.xy = rho.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type)
    rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type)
  }
  
  # p-value definition depends on the sign of the empirical correlation
  if (rho.emp>0) {
    p.perm.xy = sum(rho.null.xy>rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx>rho.emp)/nperm
  } else { 
    p.perm.xy = sum(rho.null.xy<rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx<rho.emp)/nperm
  } 
  
  # return average p-value
  return((p.perm.xy+p.perm.yx)/2)
  
}





