import pandas as pd
import nibabel as nib
import os
subdf=pd.read_csv('subid2merge.csv')
sublist=[str(n).zfill(12) for n in list(subdf['subid'])]
data_path=r'/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_PPI_1st_level'
file_name='con_0001.nii'
folder_name=['PPI_VS1_L_big-no','PPI_VS1_R_big-no']
for folder_i in folder_name:
	fu2_con_files = [os.path.exists(os.path.join(data_path,n,folder_i,file_name)) for n in sublist]
	subdf['confile']=fu2_con_files
	subdf.to_csv('FU2_test_file_exist_%s.csv' %folder_i)
	all_files = [os.path.join(data_path,n,folder_i,file_name) for n in sublist]
	concat_img=nib.funcs.concat_images(all_files)
	concat_img.to_filename('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME/merge_con_files/FU2_all_con1_%s.nii' %folder_i)

