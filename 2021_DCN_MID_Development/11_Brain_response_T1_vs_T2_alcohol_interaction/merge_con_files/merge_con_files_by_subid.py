import pandas as pd
import nibabel as nib
import os
subdf=pd.read_csv('subid2merge.csv')
sublist=[str(n).zfill(12) for n in list(subdf['subid'])]
data_path=r'/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_1st_level'
file_name='con_0002.nii'
fu2_con_files = [os.path.exists(os.path.join(data_path,n,file_name)) for n in sublist]
subdf['confile']=fu2_con_files
subdf.to_csv('BSL_test_file_exist.csv')
all_files = [os.path.join(data_path,n,file_name) for n in sublist]
concat_img=nib.funcs.concat_images(all_files)
concat_img.to_filename('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/activation_alc_LME/merge_con_files/BSL_part_con2.nii')

subdf=pd.read_csv('subid2merge.csv')
sublist=[str(n).zfill(12) for n in list(subdf['subid'])]
data_path=r'/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_1st_level'
file_name='con_0002.nii'
fu2_con_files = [os.path.exists(os.path.join(data_path,n,file_name)) for n in sublist]
subdf['confile']=fu2_con_files
subdf.to_csv('FU2_test_file_exist.csv')
all_files = [os.path.join(data_path,n,file_name) for n in sublist]
concat_img=nib.funcs.concat_images(all_files)
concat_img.to_filename('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/activation_alc_LME/merge_con_files/FU2_part_con2.nii')
