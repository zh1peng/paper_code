import numpy as np
import nibabel as nib
import os
import glob
import pandas as pd


subdf=pd.read_csv('subid2merge.csv')
sublist=[''.join(['sub-',str(n).zfill(12),'_ses-BSL']) for n in list(subdf['subid'])]
data_path=r'/gpfs2/scratch/jottinog/Freesurfer/BSL'
smoothname='fwhm25'
for hemi in ['rh','lh']:
	for surface in ['thickness']:
		file_name='.'.join([hemi,surface,smoothname,'fsaverage.mgh'])
		all_files = [os.path.exists(os.path.join(data_path,n,'surf',file_name)) for n in sublist]
		outname=".".join([hemi,'all',surface,smoothname,'mgh'])
		subdf[outname]=all_files
subdf.to_csv('BSL_exist_file.csv') 
