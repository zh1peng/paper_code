source('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/vx_LME_functions_interaction.R')

ncore2use=8
# 1 read behav data
behave_csv='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/audit_brain_cov_with_alc_group.csv'
numeric_col=c('age','mod_pds','mean_FD','age_diff')
factor_col=c('subcode','handedness','sex','all_site','time','fu2_alc')
behav_data=vx_lme_read_behav_data(behave_csv,numeric_col, factor_col)

# 2 read brain data
bsl_merged_data='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/merge_con_files/BSL_all_con1_PPI_VS1_R_big-no.nii'
fu2_merged_data='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/merge_con_files/FU2_all_con1_PPI_VS1_R_big-no.nii'

bsl_brain_data_list=vx_lme_read_merged_data(bsl_merged_data)
fu2_brain_data_list=vx_lme_read_merged_data(fu2_merged_data)

brain_data_list=list(brain_data=cbind(bsl_brain_data_list$brain_data, fu2_brain_data_list$brain_data),
		dim_info=bsl_brain_data_list$dim_info)

brain_data=brain_data_list$brain_data

# 3 true model
#formula2use='img_vx~time*fu2_alc+mod_pds+mean_FD+handedness+sex+(1|all_site/subcode)'
ptm=proc.time()
formula2use='img_vx~time*fu2_alc+mod_pds+mean_FD+handedness+sex+all_site+age_diff+(1|subcode)'
var2extract1='time2'
var2extract2='fu2_alc1'
lme_results=vx_lme_parallel_fun(brain_data, behav_data, formula2use, var2extract1, var2extract2,core_n=ncore2use)

result_path='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_alc_new_subid_balanced/PPI_R_VS'
result_filename='R_VS'
spm_path='/gpfs1/home/z/c/zcao4/matlab_tools/spm12'
#vx_lme_save_results <- function(result_list, output_path,dim_info,filename_info, spm_path)
vx_lme_save_results(lme_results, result_path,brain_data_list$dim_info,result_filename, spm_path)
print(proc.time()-ptm)







