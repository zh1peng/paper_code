
library(doParallel)
library(AnalyzeFMRI)
#library(matlabr)
#library(fslr)
library(stringr)
library(dplyr)
library(forcats)
#================================== funciton for single voxel==================================

vx_single_fun <- function (vx_i,brain_data_vx1,brain_data_vx2){
  nan_sum=sum(is.nan(brain_data_vx2))
  zero_sum=sum(brain_data_vx2==0,na.rm = TRUE)
  sub_n=length(brain_data_vx2)
  nonnan_sum=sub_n-nan_sum
  if (nan_sum<(0.5*sub_n) & zero_sum<(0.5*nonnan_sum)){ #if <half of voxels in second group are neither nan nor zero, run the following:
    g1_df=data.frame(vx_data=as.numeric(brain_data_vx1)) %>% mutate(group=1)
    g2_df=data.frame(vx_data=as.numeric(brain_data_vx2)) %>% mutate(group=2)
    tmp_data=rbind(g1_df,g2_df)
    tmp_results=var.test(vx_data~group,data=tmp_data)
    f_val=as.numeric(tmp_results$statistic)
    p_val=as.numeric(tmp_results$p.value)
    
    rm(tmp_results,tmp_data)
  }
  else{
    f_val=NaN
    p_val=NaN}
  return(cbind(f_val,p_val))
}


#====================================== Parallel implementation========================================
vx_parallel_fun <- function (brain_data1, brain_data2, core_n){
  registerDoParallel(cores=core_n)
  ptm=proc.time()
  tmp_results=foreach (vx_i = 1:dim(brain_data1)[1]) %dopar% vx_single_fun(vx_i,brain_data1[vx_i,],brain_data2[vx_i,])
  proc.time()-ptm
  f_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'f_val'])))
  p_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'p_val'])))
  result_list=list('f_vals'=f_vals,
                   'p_vals'=p_vals)
  return(result_list)
}

#==========================read data====================================
vx_read_merged_data <- function(data2read){
  # Read merged data (4D: [x,y,z]*subj) covert to voxel*subj format
  # Output:
  #   merged brain_data
  #   dim_info
  raw_brain=f.read.nifti.volume(data2read)
  n_sub=dim(raw_brain)[4]
  dim_info=dim(raw_brain)[1:3]
  brain_data=array(raw_brain,dim=c(prod(dim_info),n_sub))
  brain_data_list=list(brain_data=brain_data,
                       dim_info=dim_info)
  return(brain_data_list)
}

#======================================Save results==================================================
vx_save_results <- function(result_list, output_path,dim_info,filename_info){
  dir.create(output_path, showWarnings = FALSE)
  setwd(output_path)  
  f_write=array(result_list$f_vals,dim_info)
  p_write=array(result_list$p_vals,dim_info)
  fdrp_write=array(p.adjust(result_list$p_vals,method = 'fdr'),dim_info)
  f.write.nifti(f_write,sprintf('f_vals_%s.nii',filename_info),nii=TRUE)
  f.write.nifti(p_write,sprintf('p_vals_%s.nii',filename_info),nii=TRUE)
  f.write.nifti(fdrp_write,sprintf('fdrp_vals_%s.nii',filename_info),nii=TRUE)
}


measures=c('activiation','PPI_VS1_L_big-no','PPI_VS1_R_big-no')
title2show=c('Brain response', 'Left VS functional connectivity', 'Right VS functional connectivity')
plist=list()
for (measure2test in measures){
brain_list1=vx_read_merged_data(sprintf('/Users/zhipeng/Desktop/revision_test/BSL_g1_%s.nii',measure2test))
brain_list2=vx_read_merged_data(sprintf('/Users/zhipeng/Desktop/revision_test/BSL_g2_%s.nii',measure2test))
brain_data1=brain_list1$brain_data
brain_data2=brain_list2$brain_data
vx_results=vx_parallel_fun(brain_data1,brain_data2,10)
result_path='/Users/zhipeng/Desktop/revision_test'
result_filename=measure2test
vx_save_results(vx_results, result_path,brain_list1$dim_info,result_filename)
}


measures=c('activiation','PPI_VS1_L_big-no','PPI_VS1_R_big-no')
title2show=c('Brain response', 'Left VS functional connectivity', 'Right VS functional connectivity')
plist=list()
for ( m_i in c(1:3)){
  measure2test=measures[m_i]
  brain_list1=vx_read_merged_data(sprintf('/Users/zhipeng/Desktop/revision_test/BSL_g1_%s.nii',measure2test))
  brain_list2=vx_read_merged_data(sprintf('/Users/zhipeng/Desktop/revision_test/BSL_g2_%s.nii',measure2test))
  brain_data1=brain_list1$brain_data
  brain_data2=brain_list2$brain_data
  vx_results=vx_parallel_fun(brain_data1,brain_data2,10)
  df2plot=data.frame(f_vals=unlist(vx_results$f_vals),
                     p_vals=unlist(vx_results$p_vals))
  p=ggplot(df2plot, aes(x=f_vals))+
    geom_histogram()+
    xlab('F value')+
    ggtitle(title2show[m_i])
  plist[[m_i]]=p
}
main=grid.arrange(grobs=plist, nrow = 1)
ggsave('/Users/zhipeng/Google Drive/post-doc/MID_network_FU2/process_code/intersect_cov/AUDIT analysis new subid/revision_analysis_drop_out_test/F_dist.tiff',
       main,width = 10,height = 4)




