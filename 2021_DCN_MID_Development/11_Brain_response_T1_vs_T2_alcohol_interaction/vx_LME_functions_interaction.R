# handy function used for vx_LME analysis
# timeline log:
# 2020/08/13
# zhipeng


library(lmerTest)
library(doParallel)
library(AnalyzeFMRI)
library(matlabr)
#library(fslr)
library(stringr)
library(dplyr)
library(forcats)
#================================== LME funciton for single voxel==================================

vx_lme_single_fun <- function (vx_i,brain_data_vx,behav_data,formula2use,var2extract1,var2extract2){
# vx_lme_single_fun:
# apply LME model on the voxel/vertix/edgewise imaging data
# allows parallel computation
# vx_i: the ith voxel
# brain_data_vx: the voxel
# behav_data: cov inlcuded in the LME model
# formula2use: LME formula2use
# var2extract1: var1
# var2extract2: var2
# the interaction of var1 and var2 will be extracted as well

  tmp_data=behav_data
  nan_sum=sum(is.nan(brain_data_vx))
  zero_sum=sum(brain_data_vx==0,na.rm = TRUE)
  sub_n=length(brain_data_vx)
 nonnan_sum=sub_n-nan_sum
  if (nan_sum<(0.5*sub_n) & zero_sum<(0.5*nonnan_sum)){ #if <half of voxels are neither nan nor zero, run the following:
    tmp_data['img_vx']=as.numeric(brain_data_vx)
    tmp_lme_model=lmer(formula(formula2use),data=tmp_data, na.action=na.omit,
                       control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    t1_val=summary(tmp_lme_model)$coefficients[var2extract1,'t value']
    p1_val=summary(tmp_lme_model)$coefficients[var2extract1,'Pr(>|t|)']
    t2_val=summary(tmp_lme_model)$coefficients[var2extract2,'t value']
    p2_val=summary(tmp_lme_model)$coefficients[var2extract2,'Pr(>|t|)']
    int_t_val=summary(tmp_lme_model)$coefficients[sprintf('%s:%s',var2extract1,var2extract2),'t value']
    int_p_val=summary(tmp_lme_model)$coefficients[sprintf('%s:%s',var2extract1,var2extract2),'Pr(>|t|)']
    rm(tmp_lme_model,tmp_data)
  }
  else{
    t1_val=NaN
    p1_val=NaN
    t2_val=NaN
    p2_val=NaN
    int_t_val=NaN
    int_p_val=NaN}
  return(cbind(t1_val,p1_val,t2_val,p2_val,int_t_val,int_p_val))
}


#======================================LME Parallel implementation========================================
vx_lme_parallel_fun <- function (brain_data, behav_data, formula2use, var2extract1,var2extract2, core_n){
# Input:
# brain_data: voxels * subj
# behav_data: behav variables and covirates
# formular2use: LME fomula (img_vx is brain voxel)
#               img_vx~S5+sex_2+FU2_age+handedness+C1+C2+C3+C4+(1|all_site)
# var2extract: what predictor you want to get the beta/t/p values
# core_n: how many cores you want to use in the parallel processing
# 
# Output:
# result_list: beta/t/p values corresponding to var2extract

registerDoParallel(cores=core_n)
ptm=proc.time()
tmp_results=foreach (vx_i = 1:dim(brain_data)[1]) %dopar% vx_lme_single_fun(vx_i,brain_data[vx_i,],behav_data,formula2use,var2extract1,var2extract2)
proc.time()-ptm
t1_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'t1_val'])))
p1_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'p1_val'])))
t2_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'t2_val'])))
p2_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'p2_val'])))
t1_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'t1_val'])))
p1_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'p1_val'])))
int_t_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'int_t_val'])))
int_p_vals=unlist(lapply(tmp_results, function(i) as.numeric(i[,'int_p_val'])))

result_list=list('t1_vals'=t1_vals,
               'p1_vals'=p1_vals,
               't2_vals'=t2_vals,
               'p2_vals'=p2_vals,
               'int_t_vals'=int_t_vals,
               'int_p_vals'=int_p_vals)
return(result_list)
}

#======================================Prepare behav_data===============================================  
vx_lme_read_behav_data <- function (behav_csv,numeric_col,factor_col){
# behav_csv: behavioural cvs to read
# numeric_col: e.g., age S5
# factor_col: e.g., site, hadedness
  
  all_vars=read.csv(behav_csv)
  for(col in factor_col){
    all_vars[,col]=as.factor(as.character(all_vars[,col]))
  }
  for (col in numeric_col){
    all_vars[,col]=scale(as.numeric(as.character(all_vars[,col])),center=TRUE, scale=TRUE)
  }
  col2use=append(numeric_col,factor_col)
  behav_data=all_vars[,col2use]
  return(behav_data)
}

#======================================Prepare brain_data===============================================
vx_lme_read_subject_nii <- function(nii2read){
# Read subject-level data and merge them as voxel*subj format
# Only works with the strucutre like: data_path/subfolder[i]/nii_name
# Input:
#     nii2read: a char vector of data_path/subfolder[i]/nii_name
# Output:
#     merged brain_data
#     dim_info

  brain_data=numeric()
  for (subi in 1:length(nii2read)){
    sub_brain=f.read.nifti.volume(nii2read[subi])
    dim_info=dim(sub_brain)[1:3]
    brain_data=cbind(brain_data,array(sub_brain,dim=c(prod(dim_info))))
  }
  
  brain_data_list=list(brain_data=brain_data,
                       dim_info=dim_info) 
  return(brain_data_list)
} 

vx_lme_read_merged_data <- function(data2read){
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

matlab_tfce<-function(nii_file,spm_path,filename_info){
matlab_code=c(sprintf("addpath('%s')",spm_path),
     sprintf("vol_info=spm_vol('%s')",nii_file),
  "img=spm_read_vols(vol_info)",
sprintf("mask_info=spm_vol('%s/matched_GM_mask.nii')",spm_path),
"mask_img=spm_read_vols(mask_info)",
"tmp_img=img.*mask_img",
"tmp_img(tmp_img<0)=0",
  "[pos_tfced] = matlab_tfce_transform(tmp_img,2,0.5,26,0.1)",
sprintf("mask_info.fname='t_vals_%s_tfced_pos.nii'",filename_info),
  "spm_write_vol(mask_info,pos_tfced)",
"tmp_img=img.*mask_img",
"tmp_img(tmp_img>0)=0",
  "[neg_tfced] = matlab_tfce_transform(abs(tmp_img),2,0.5,26,0.1)",  
  sprintf("mask_info.fname='t_vals_%s_tfced_neg.nii'",filename_info),
  "spm_write_vol(mask_info,neg_tfced)")
res=run_matlab_code(matlab_code)
}



#======================================Save results==================================================
vx_lme_save_results <- function(result_list, output_path,dim_info,filename_info, spm_path){
# write beta/t/p values into nifiti
# why use matlab?
# spm/xjview cannot read the output from R. So I read-->write the files again using SPM and match to a spm templeate file
# 
# Input:
#   result_list: generate from vl_lme_parallel_fun
#   output_path: where to save the files
#   dim_info: [x,y,z] dimension info to covert data back
#   filename_info: what info to include in the filename
#   spm_path: your spm_path should include this code matlab_tfce_transform.m !!!!
dir.create(output_path, showWarnings = FALSE)
setwd(output_path)  
  
t1_write=array(result_list$t1_vals,dim_info)
p1_write=array(result_list$p1_vals,dim_info)
t2_write=array(result_list$t2_vals,dim_info)
p2_write=array(result_list$p2_vals,dim_info)
int_t_write=array(result_list$int_t_vals,dim_info)
int_p_write=array(result_list$int_p_vals,dim_info)


f.write.nifti(t1_write,sprintf('t1_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(p1_write,sprintf('p1_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(t2_write,sprintf('t2_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(p2_write,sprintf('p2_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(int_t_write,sprintf('int_t_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(int_p_write,sprintf('int_p_vals_%s.nii',filename_info),nii=TRUE)
}





#======================================Save null results==================================================
vx_lme_save_null_results <- function(result_list, output_path,dim_info,filename_info, spm_path){
# write beta/t/p values into nifiti
# why use matlab?
# spm/xjview cannot read the output from R. So I read-->write the files again using SPM and match to a spm templeate file
# 
# Input:
#   result_list: generate from vl_lme_parallel_fun
#   output_path: where to save the files
#   dim_info: [x,y,z] dimension info to covert data back
#   filename_info: what info to include in the filename
#   spm_path: your spm_path should include this code matlab_tfce_transform.m !!!!
dir.create(output_path, showWarnings = FALSE)
setwd(output_path)  
t1_write=array(result_list$t1_vals,dim_info)
t2_write=array(result_list$t2_vals,dim_info)
int_t_write=array(result_list$int_t_vals,dim_info)
f.write.nifti(t1_write,sprintf('t1_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(t2_write,sprintf('t2_vals_%s.nii',filename_info),nii=TRUE)
f.write.nifti(int_t_write,sprintf('int_t_vals_%s.nii',filename_info),nii=TRUE)
}




#======================================extract permutation==================================================
vx_lme_extract_perm_results <- function(perm_result_path){
# read permutation results
# Input:  perm_result_path
#         file_pattern '^t_vals_shuffle*tfced\\.txt$'
# Output: sorted TFCE max values
setwd(perm_result_path)
tfce2read=list.files(path = perm_result_path, pattern = 'tfced_pos\\.nii$',no.. = TRUE)       
tmp_max=lapply(tfce2read, function(x) max(f.read.nifti.volume(x)))
perm_tfce_max=sort(unlist(tmp_max)) # sorted
write.csv(perm_tfce_max,'perm_tfce_pos_max.csv')

tfce2read=list.files(path = perm_result_path, pattern = 'tfced_neg\\.nii$',no.. = TRUE)       
tmp_max=lapply(tfce2read, function(x) max(f.read.nifti.volume(x)))
perm_tfce_max=sort(unlist(tmp_max)) # sorted
write.csv(perm_tfce_max,'perm_tfce_neg_max.csv')


tvals2read=list.files(path = perm_result_path, pattern = '\\d.nii$',no.. = TRUE)       
tmp_max=lapply(tvals2read, function(x) max(abs(f.read.nifti.volume(x)),na.rm=TRUE))
perm_tvals_max=sort(unlist(tmp_max)) # sorted
write.csv(perm_tvals_max,'perm_tvals_max.csv')
return(perm_tfce_max)
}


#======================================calculate P values==================================================
vx_lme_calculate_fwep <- function(perm_tfce_max,result_path,result_filename){
  setwd(result_path) 
  p_array=(length(perm_tfce_max):1)/length(perm_tfce_max)
  result_file=sprintf('t_vals_%s_tfced.nii',result_filename)
  raw_tfce=f.read.nifti.volume(result_file)
  tmp_p=sapply(raw_tfce,function(x) p_array[which.min(abs(perm_tfce_max-x))]) # perm_tfce_max is vecotr of max tfce values from permutation
  corrected_p=array(tmp_p,dim(raw_tfce)[1:3])
  f.write.nifti(corrected_p,sprintf('p_vals_fwe_%s.nii',result_filename),nii=TRUE)
}


vx_lme_calculate_fdr <- function(result_path, result_filename){
t_corrected=t2write[t2write<t_fdr_df_max]
df_max=round(max(result_list$df_vals,na.rm = TRUE))
t_fdr_df_max=Threshold.FDR(t2write,0.05,df1=df_max)
return(t_fdr_df_max)
}





#======================================shuffle time and alc label in cov==================================================
shuffle_time_alc <- function(cov){
  bsl_df=cov%>% filter(time==1)
   tmp_fu2_alc=sample(bsl_df$fu2_alc) # alc is same for bsl and fu2, so it is shuffled between subjects
     cov$shuffled_fu2_alc=fct_c(tmp_fu2_alc,tmp_fu2_alc)
  # for the subjects that are shown in both BSL and FU2
  # shuffle within subjects to prevent same subject in same wave
  cov$shuffled_time=cov$time # there seems no need to shuffle time
  return(cov)
}
