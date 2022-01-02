#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:18:40 2020

@author: zhipeng
"""


import glob
import os
import nibabel as nib
output = r'/gpfs2/scratch/zcao4/IMAGEN_MID_data/MID_FU2_EPI'
search_test = r'/gpfs2/scratch/acjulian/IMAGEN_files_still_need_BIDSifying/Task_FU2_Funcs/*/EPI_mid/*sw*.nii.gz'
#len(list(glob.iglob(search_test, recursive=True)))
if not os.path.exists(output):
    os.makedirs(output)
count = 0
wrong_subid = []
all_subid = []
for filen in glob.iglob(search_test, recursive=True):
    subid = os.path.basename(os.path.dirname(os.path.dirname(filen)))
    subpath = os.path.join(output, subid)
    if not os.path.exists(subpath):
        os.makedirs(subpath)
    count += 1
    print('=================' + str(count) + '==================')
    try:  # split 4d to 3d to directory
        img = nib.load(filen)
        imgs = nib.four_to_three(img)
        froot, ext = os.path.splitext(filen)
        froot, ext = os.path.splitext(froot)
        _, fname = os.path.split(froot)
        froot = os.path.join(subpath, fname)
        for i, img3d in enumerate(imgs):
            fname3d = '%s_%05d.nii' % (froot, i + 1)
            nib.save(img3d, fname3d)
        all_subid.append(subid)
    except:
        print(subid + 'seems not correct')
        wrong_subid.append(subid)

with open(os.path.join(output, 'subid.txt'), 'w') as file:
    for line in all_subid:
        file.write(line)
        file.write('\n')
with open(os.path.join(output, 'log.txt'), 'w') as file:
    for line in wrong_subid:
        file.write(line)
        file.write('\n')