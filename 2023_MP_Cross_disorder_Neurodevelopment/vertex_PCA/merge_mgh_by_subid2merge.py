# -*- coding: utf-8 -*-
"""
merge mgh by subid2merge

"""
import numpy as np
import nibabel as nib
import os
import glob
import pandas as pd
def loadmgh(imagename):
	if os.path.exists(imagename): # check if file exists
		img = nib.freesurfer.mghformat.load(imagename)
		img_data = img.get_data()
	else:
		print("Cannot find input image: %s" % imagename)
		exit()
	return (img,img_data)

def savemgh(imgdata, img, index, imagename):
	outdata = imgdata.astype(np.float32, order = "C")
	if imgdata.ndim == 2:
		imgout = np.zeros((img.shape[0],img.shape[1],img.shape[2],outdata.shape[1]))
	elif imgdata.ndim == 1:
		imgout = np.zeros((img.shape[0],img.shape[1],img.shape[2]))
	else:
		print('error')
	imgout[index]=outdata
	nib.save(nib.freesurfer.mghformat.MGHImage(imgout.astype(np.float32, order = "C"),img.affine),imagename)


# test if mgh files for subid2merge exist
""" subdf=pd.read_csv('subid2merge.csv')
sublist=[''.join(['sub-',str(n).zfill(12),'_ses-FU2']) for n in list(subdf['subid'])]
data_path=r'/gpfs2/scratch/jottinog/Freesurfer/FU2'
smoothname='fwhm25'
for hemi in ['rh','lh']:
	for surface in ['thickness']:
		file_name='.'.join([hemi,surface,smoothname,'fsaverage.mgh'])
		all_files = [os.path.exists(os.path.join(data_path,n,'surf',file_name)) for n in sublist]
		outname=".".join([hemi,'all',surface,smoothname,'mgh'])
		subdf[outname]=all_files
subdf.to_csv('FU2_exist_file.csv')  """



data_path=r'/gpfs2/scratch/jottinog/Freesurfer/BSL'
subdf=pd.read_csv('subid2merge.csv')
sublist=[''.join(['sub-',str(n).zfill(12),'_ses-BSL']) for n in list(subdf['subid'])]
for smoothname in ['fwhm0','fwhm5','fwhm10','fwhm15','fwhm20','fwhm25']:
    for hemi in ['rh','lh']:
        for surface in ['thickness']:
            file_name='.'.join([hemi,surface,smoothname,'fsaverage.mgh'])
            all_files = [os.path.join(data_path,n,'surf',file_name) for n in sublist]
            outname=".".join([hemi,'all',surface,smoothname,'mgh'])
            numMerge=len(all_files)
            img, img_data = loadmgh(all_files[0])
            mask_index = np.zeros((img_data.shape[0],img_data.shape[1],img_data.shape[2]))
            mask_index = (mask_index == 0)
            img_data_trunc = img_data[mask_index].astype(np.float32)
            for i in range(numMerge):
                    print("merging image %s" % all_files[i])
                    if i > 0:
                        _, tempimgdata = loadmgh(all_files[i])
                        tempimgdata=tempimgdata[mask_index].astype(np.float32)
                        img_data_trunc = np.column_stack((img_data_trunc,tempimgdata))
                # remove nan if any
            img_data_trunc[np.isnan(img_data_trunc)] = 0
            savemgh(img_data_trunc.astype(np.float32), img, mask_index, outname)

 
