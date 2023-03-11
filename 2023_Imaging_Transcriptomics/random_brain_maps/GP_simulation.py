#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:26:00 2022

@author: zhipeng
"""
import os
import numpy as np
from scipy import ndimage
from netneurotools.freesurfer import _decode_list
import pandas as pd
# prob need to save the file as we need to look at different atlas


dir2save='/Users/zhipeng/Desktop/simpath'
n_sim=1000

hemi='lh'
anot_file='/Applications/freesurfer/subjects/fsaverage5/label/{}.aparc.annot'.format(hemi)
drop=['unknown', 'corpuscallosum']
labels, ctab, names = nib.freesurfer.read_annot(anot_file)
names=_decode_list(names)
affine = np.eye(4) * 2
affine[:, -1] = [-90, -90, -72, 1]

for alpha in (0.1,0.2,0.3):
    print(alpha)
    alpha_folder=os.path.join(dir2save,'alpha_'+str(alpha))
    os.mkdir(alpha_folder)

    for seed_i in range(n_sim):
        print(seed_i)
        gfield = gaussian_random_field(91, 109, 91, alpha=alpha, seed=seed_i)
        fn = make_tmpname(suffix='.nii.gz')
        nib.save(nib.nifti1.Nifti1Image(gfield, affine), fn)
        outname = os.path.join(alpha_folder,hemi+'_simbrain'+str(seed_i).zfill(3)+'.mgh')
        run(VOL2SURF.format(fn, outname, hemi), quiet=True)
        os.remove(fn)
    
    tmp=list()
    colname=list()
    for seed_i in range(n_sim):
        mgh_file= os.path.join(alpha_folder,hemi+'_simbrain'+str(seed_i).zfill(3)+'.mgh')         
        img = nib.load(mgh_file).dataobj
        data = ndimage.mean(np.squeeze(img), labels, np.unique(labels))
        drop = np.intersect1d(names, drop)
        data = np.delete(data, [names.index(f) for f in drop])            
        tmp.append(data)
        colname.append('simbrain'+str(seed_i).zfill(3))    
       
    df=pd.DataFrame(np.array(tmp).T, index =['L_'+f for f in names if f not in drop],
                                                  columns =colname)
    df.to_csv(os.path.join(alpha_folder,'alpha'+str(alpha)+'_'+hemi+'_simbrain_DK.csv'),  index=True)









