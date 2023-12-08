# all functions included in the analysis
#===============================================
# 1. preprocess_df
# preprocess the dataframe including:
# a. set all data to the right 
# b. rescale ICV
# c. find the largest site and set it to the reference
# d. set outlier to 3*IQR
preprocess_df <- function(df,Group=F){
  #if Group==T set Group variable as factor.
  
  #set data format
  table(apply(df, 2, class)) ## all columns are of class character
  
  ##===========================Set numeric cols==============================================
  df2 <- df
  for(col in colnames(df2)[grep("^L_|^R_|Age|ICV|mod_PDS", colnames(df))] ){
    df2[,col] <- as.character(df[,col])
    df2[,col] <- as.numeric(df2[,col])
  }
  ## checked, values don't change after setting to numeric. FALSE means no changes
  TRUE %in% apply(df==df2, 2, function(x){FALSE %in% x}) 
  df <- df2
  rm(col, df2)
  
  df$ICV=df$ICV/1000000 #rescale ICV
  
  df$Sex=as.factor(as.character(df$Sex))
  df$Sex=relevel(df$Sex,'0')
  
  df$Site=as.factor(as.character(df$Site))
  fac_num=df %>% 
    group_by(Site) %>%
    dplyr::summarise(n_rows = length(Site))%>%
    arrange(n_rows)%>% 
    mutate_if(is.factor, as.character) # as.charater(factor) will give only numbers!!! as.character(fac_num[dim(fac_num)[1],1])
  # set the largest group as reference
  df$Site=relevel(df$Site,as.character(fac_num[dim(fac_num)[1],1]))
  
  if (Group){
    df$Group=as.factor(as.character(df$Group))
    df$Group=relevel(df$Group,'0')
  }
  
  ##=============================set 0s as NA============================================
  for(i in colnames(df)[grep("L_|R_|ICV", colnames(df))] ){
    if(length(which(df[,i]==0)) !=0){
      print(paste0("remove ",length(which(df[,i]==0))," NULLs for ",i))
      df[which(df[,i]==0),i] <- NA
      df[,i]=as.numeric(df[,i])
    }
  }
  df=df%>%filter(complete.cases(.))
  
  #outlier clipper
  for(region in colnames(df)[grep("L_|R_", colnames(df))]){
    
    #Outliers were identified as above or below three times k=3
    # the interquartile range
    # https://statisticsbyjim.com/basics/outliers/
    
    k=3
    Q1=as.numeric(quantile(df[,region],na.rm = TRUE)[2])
    Q3=as.numeric(quantile(df[,region],na.rm = TRUE)[4])
    IQR_k=k*(Q3-Q1)
    low_fence=Q1-IQR_k
    up_fence=Q3+IQR_k
    
    #anoter approach is 3std
    #low_fence=mean(df[,region],na.rm = TRUE)-3*as.numeric(sd(df[,region],na.rm = TRUE))
    #up_fence=mean(df[,region],na.rm = TRUE)+3*as.numeric(sd(df[,region],na.rm = TRUE))
    
    ## clip
    df[which(df[,region]>up_fence),region] <- up_fence
    df[which(df[,region]<low_fence),region] <- low_fence
  }
  return(df)
}


#===============================================
# 2. resid_CT
# do resid for CT
# type=response raw residuals (no scaling applied)
resid_CT <- function(df,
                     if_pds=F,
                     if_meanCT=F){
  # resid sex/age/ICV/site
  # if_pds=T, also include mod_PDS in the formular
  # if_meanCT=T, also include overall mean CT in the formular
  lme_models=list()
  resid_data=list()
  df$meanCT=rowMeans(df[grep("^L_|^R_", colnames(df))])
  base_fomular='outcome~Sex + Age + ICV'
  if (if_pds){
    base_fomular=paste(base_fomular,'mod_PDS',sep='+')
  }
  
  if (if_meanCT){
    base_fomular=paste(base_fomular,'meanCT',sep='+')
  }
  
  for (region in colnames(df)[grep("^L_|^R_", colnames(df))] ){
    df$outcome <- df[,region]
    lme_models[[region]] <- lme(as.formula(base_fomular), random=~1|Site, data = df,
                                method='REML',control = lmeControl(opt = "optim"))  
    resid_data[[region]] <- resid(lme_models[[region]], type='response') #type=response raw residuals
  }
  resid_df=data.frame(resid_data)%>% rename_all( ~ paste0("resid_", .x))
  return(resid_df)
}



resid_CT_by_formula <- function(df, f2use='outcome~Sex+Age+ICV'){
  # site is added as random by default
  # outcome represents the L/R regions
  lme_models=list()
  resid_data=list()
  for (region in colnames(df)[grep("^L_|^R_", colnames(df))] ){
    df$outcome <- df[,region]
    lme_models[[region]] <- lme(as.formula(f2use), random=~1|Site, data = df,
                                method='REML',control = lmeControl(opt = "optim"))  
    resid_data[[region]] <- resid(lme_models[[region]], type='response') #type=response raw residuals
  }
  resid_df=data.frame(resid_data)%>% rename_all( ~ paste0("resid_", .x))
  return(resid_df)
}










#===============================================
# 3. do_pca_df
do_pca_df <- function(df,
                      if_resid=T,...){
  # PCA pipeline
  if (if_resid){
    resid_df=resid_CT(df,...)
    pca_data=as.matrix(resid_df)}
  else {pca_data=as.matrix(df[,grep("^L_|^R_", colnames(df))])}
  raw_data=as.matrix(df[,grep("^L_|^R_", colnames(df))])
  pca<-prcomp(pca_data, center = T, scale. = T)
  #adjust sign according to the diagnal: make all diag(U)>0
  rotation.orig=pca$rotation
  pca$rotation=pca$rotation%*%diag(as.numeric(diag(rotation.orig)>0)*2-1) # keep pos diag eig but flip neg diag eig col with -1
  pca$x=pca$x%*%diag(as.numeric(diag(rotation.orig)>0)*2-1) 
  pca$stand_loadings=stand_loadings(pca, pca_data) # function from syndRomics
  pca$pca_data=pca_data
  pca$raw_data=raw_data
  pca$raw_df=df
  pca=replace_loading_names(pca) #rename long names to short names for barplot
  return(pca)
}

# 4. do_pca for a dataset
do_pca <- function(data_path='F:/Google Drive/post-doc/p_factor',
                   data_name='IMAGEN_FU2',...){
  # the pca pipeline
  # if_resid=T do resid before pca
    # if_pds if need to add pds to resid
    # if_meanCT if need to  add meanCT to resid
  
  # PCA pipeline
  data_file=sprintf('%s/%s/%s_CT_ready.csv',data_path,data_name,data_name)
  #read data
  df=read.csv(data_file)
  df=preprocess_df(df,Group = F)
  pca=do_pca_df(df,...)
  return(pca)
}

#===============================================
# 5 repalce loading names with abbreviations for display purpose
# need to specify the dir to dk_abbrevation file
replace_loading_names <- function(pca.obj,
                                  keyfile='F:/Google Drive/post-doc/p_factor/clean_codes/dk_abbrevation.csv'){
  
  
  key_df=read.csv(keyfile)
  patterns=as.character(key_df$region)
  replacement=as.character(key_df$abr)
  
  load_df=pca.obj$stand_loadings
  row.name=rownames(load_df)
  row.name=stringr::str_replace_all(row.name,setNames(replacement,patterns))
  row.name=stringr::str_replace_all(row.name,'_thickavg|resid_','')
  row.name=stringr::str_replace_all(row.name,'_','.')
  renamed_loading=load_df
  rownames(renamed_loading) <- row.name
  pca.obj$renamed_loading=renamed_loading
  return (pca.obj)
}


#===============================================
# 6. Remove pc1 from the data
# data is reconstructed using PCs without PC1
# results are stored in pca.obj$pc1_rm_data
remove_pc1 <- function(pca.obj){
  nComp=ncol(pca.obj$x)
  pca.obj$pc1_rm_data=t(t(pca.obj$x[,2:nComp] %*% t(pca.obj$rotation[,2:nComp])) * pca.obj$scale + pca.obj$center)
  pca.obj$pc1_rm_data=as.data.frame(pca.obj$pc1_rm_data) %>% rename_with(~paste0('pc1_rm_',.))
  return(pca.obj)
}

keep_pc1 <- function(pca.obj){
  pca.obj$pc1_only_data=t(t(pca.obj$x[,1] %*% t(pca.obj$rotation[,1])) * pca.obj$scale + pca.obj$center)
  pca.obj$pc1_only_data=as.data.frame(pca.obj$pc1_only_data) %>% rename_with(~paste0('pc1_only_',.))
  return(pca.obj) 
}


#===============================================
# 7. make a plot to show brain/brain/VAF
plot_pca <- function(pca.obj,title2add='title'){
  p1= VAF_plot(pca.obj, ndim =1:10)+ylab('VAF (%)')
  loadings=pca.obj$stand_loadings
  p2=loadings%>%
    select(PC1,PC2,PC3)%>%
    mutate(label=rownames(loadings),
           label=sub('.*L_','lh_',label),
           label=sub('.*R_','rh_',label),
           label=sub('_thickavg','',label))%>%
    gather(PCs,Loadings, PC1:PC3)%>%
    group_by(PCs) %>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = Loadings)) +
    geom_sf(show.legend = F) +
    facet_wrap( ~ PCs,ncol=1)+
    scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=c(-1,1),aes(title="Loadings"))+
    theme_void()
  p3=barmap_loading_rename(pca.obj, ndim=1:3, star_values = F,text_values = T)
  lay_out=rbind(c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(1,1,1,2,2,2),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3),
                c(3,3,3,3,3,3))
  main=grid.arrange(p1,p2,p3,top=title2add,layout_matrix = lay_out)
}

#===============================================
# 8. just plot loadings
plot_pca_brain <- function (pca.obj,pc2plot=1){
  loadings=pca.obj$stand_loadings
  brain_p=pca.obj$stand_loadings%>%
    select(sprintf('PC%d',pc2plot))%>%
    mutate(label=rownames(loadings),
           label=sub('.*L_','lh_',label),
           label=sub('.*R_','rh_',label),
           label=sub('_thickavg','',label))%>%
    rename('loadings'=sprintf('PC%d',pc2plot))%>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = loadings)) +
    geom_sf(show.legend = F) +
    scale_fill_gradient2(midpoint=0, low="steelblue1", mid="white",high="firebrick1", space ="Lab" ,limits=c(-1,1),aes(title="Loadings"))+
    theme_void()
  return(brain_p)
}


plot_pca_brain2 <- function (pca.obj,pc2plot=1,showlegend=F){
  loadings=pca.obj$stand_loadings
  brain_p=pca.obj$stand_loadings%>%
    select(sprintf('PC%d',pc2plot))%>%
    mutate(label=rownames(loadings),
           label=sub('.*L_','lh_',label),
           label=sub('.*R_','rh_',label),
           label=sub('_thickavg','',label))%>%
    rename('loadings'=sprintf('PC%d',pc2plot))%>%
    brain_join(dk) %>%
    reposition_brain(.~hemi+side) %>%
    ggplot(aes(fill = loadings)) +
    geom_sf(show.legend = showlegend) +
    scale_fill_gradient2(midpoint=0.335, low="steelblue1", mid="darkseagreen1",high="firebrick1", space ="Lab" ,limits=c(-0.15,0.82),aes(title="Loadings"))+
    #scale_fill_gradient( low="green",high="firebrick1", space ="Lab" ,limits=c(0,0.8),aes(title="Loadings"))+
    #viridis::scale_fill_viridis(limits = c(0, 0.8))+
    #scale_fill_gradientn(colors = viridis_pal()(9), limits=c(0, 0.8))+
    #(low="green", high="red",limits=c(0, 1))+
    theme_void()
  return(brain_p)
}
#plot_pca_brain2(pca.obj,pc2plot = pc2plot,T)




#===============================================
# 9. plot VAF for a list of pca resutls
plot_pca_list_VAF <- function(pca_list,pc2plot=8,type='bar'){
  original_VAF=list()
  for (data_i in names(pca_list)){
    pca.obj=eval(parse(text=sprintf('pca_list$%s',data_i)))
    original_VAF[[data_i]]=pca.obj$sdev^2/sum(pca.obj$sdev^2)
  }
  plot_df=data.frame(original_VAF) %>% 
    slice_head(n=pc2plot) %>% 
    mutate(PC_N=paste0('PC',c(1:pc2plot),sep='')) %>% 
    gather('Data','VAF','ABCD_T1':'ENIGMA_HC') %>% 
    mutate(Data=sub('ENIGMA_HC','ENIGMA-CTRL',Data),
           Data=sub('IMAGEN_','IMAGEN-',Data),
           Data=sub('ABCD_','ABCD-',Data))
  plot_df$PC_N=factor(plot_df$PC_N,levels = paste0('PC',c(1:pc2plot),sep=''))
  if (type=='bar'){
    p=ggplot(plot_df,aes(fill=Data,x=PC_N,y=.data$VAF*100))+
      geom_bar(position="dodge", stat="identity")+
      scale_fill_manual("legend", values = c("#f38905",
                                             "#00b9d5",
                                             "#a40747",
                                             "#7ca600",
                                             "#d0bcff",
                                             "#006846"))+
    
      xlab(NULL)+ylab('VAF (%)')+ 
      theme_minimal()+
      theme(legend.position = c(0.8,0.8), 
            legend.title = element_blank(),
            axis.text.x=element_text(size=20),
            axis.text.y=element_text(size=20),
            text=element_text(size=20)
            )}
  else if (type=='line'){
    p=ggplot(plot_df,aes(.data$PC_N,.data$VAF*100,color=.data$Data))+
      geom_point(show.legend = T)+
      geom_line(aes(group=.data$Data), show.legend = T)+
      theme_minimal()+
      xlab(NULL)+ylab('VAF (%)')+
      theme(legend.position = c(0.8,0.8), 
            legend.title = element_blank(),
            text=element_text(size=16))+
      scale_color_brewer(palette="Accent")
  }
  return(p)
}



### plot related functions
replace_loading_names <- function(pca.obj,
                                  keyfile='F:/Google Drive/post-doc/p_factor/analysis_codes/dk_abbrevation.csv'){
  #repalce loading names with abbreviations for display purpose
  
  key_df=read.csv(keyfile)
  patterns=as.character(key_df$region)
  replacement=as.character(key_df$abr)
  
  load_df=pca.obj$stand_loadings
  row.name=rownames(load_df)
  row.name=stringr::str_replace_all(row.name,setNames(replacement,patterns))
  row.name=stringr::str_replace_all(row.name,'_thickavg|resid_','')
  row.name=stringr::str_replace_all(row.name,'_','.')
  renamed_loading=load_df
  rownames(renamed_loading) <- row.name
  pca.obj$renamed_loading=renamed_loading
  return (pca.obj)
}

