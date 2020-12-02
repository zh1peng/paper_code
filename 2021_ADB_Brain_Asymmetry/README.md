

# Mapping cortical and subcortical asymmetries in substance dependence: Findings from the ENIGMA Addiction Working Group

**Some of the analysis code was taken from:**
Postema, M. C., Van Rooij, D., Anagnostou, E., Arango, C., Auzias, G., Behrmann, M., ... & Deruelle, C. (2019). Altered structural brain asymmetry in autism spectrum disorder in a study of 54 datasets. Nature communications, 10(1), 1-12. [Code Link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-13005-8/MediaObjects/41467_2019_13005_MOESM4_ESM.zip) 


I have assembled the analysis code into one file: **0_run_all_analysis**:
```R
working_dir <- 'F://Google Drive//post-doc//Laterality//Manuscript//Brain_Asymmetry_upload_ADB//code2upload//'
setwd(working_dir)

##=================================Model 1=============================================
source(paste(working_dir,"1_load_data.R", sep=""))           
source(paste(working_dir,"2_model1_analysisl.R", sep=""))        
source(paste(working_dir,"3_model1_bar_plot.R", sep=""))
source(paste(working_dir,"4_model1_cohend_map.R",sep=""))


##=================================Model 2=============================================
rm(list=setdiff(ls(), "working_dir"))
for (substance2test in c('Substance_name_alc',
                      'Substance_name_nic',
                      'Substance_name_coc',
                      'Substance_name_met',
                      'Substance_name_can')){
  source(paste(working_dir,"1_load_data.R", sep=""))
  source(paste(working_dir,"5_model2_analysis.R", sep=""))
  source(paste(working_dir,"6_model2_cohend_map.R", sep=""))
  rm(list=ls()[! ls() %in% c("working_dir","substance2test")])}


# make multiple comparision correction across substance types for model2
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"7_model2_mcc_across_substance.R", sep=""))

# comparision between substance groups for the significant regions
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"1_load_data.R", sep=""))
source(paste(working_dir,'8_model2_post_hoc_analysis.R',sep=""))


# make bar plot for model 2
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"1_load_data.R", sep=""))
source(paste(working_dir,'9_model2_bar_plot.R',sep=""))


# make violin plot for model 1 and model 2
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"1_load_data.R", sep=""))
source(paste(working_dir,"10_residulize_aidata.R", sep=""))
source(paste(working_dir,"11_residulize_LRdata.R", sep=""))
source(paste(working_dir,'12_residulized_plots.R',sep="")) # model 1 violin plot


##================================= Machine learning SVM (python) ==========================================
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"13_ML_load_data.R", sep=""))
source(paste(working_dir,"10_residulize_aidata.R", sep=""))
write.csv(resid_aidata,'residulized_aidata.csv')
## run ML model using python
library(reticulate)
py_run_file("svm_classification_code.py")

##==========================Past 30 days analysis==============================
rm(list=setdiff(ls(), "working_dir"))
source(paste(working_dir,"1_load_data.R", sep=""))
source(paste(working_dir,"14_past_30_day_alc_analysis.R", sep=""))
source(paste(working_dir,"15_past_30_day_nic_analysis.R", sep=""))

```


**The outputs can be found in [output_files](https://github.com/zh1peng/paper_code/tree/main/2021_ADB_Brain_Asymmetry/output_files)**

