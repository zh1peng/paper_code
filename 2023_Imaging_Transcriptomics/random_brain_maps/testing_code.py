#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 10:26:00 2022

@author: zhipeng
"""

# prob need to save the file as we need to look at different atlas
affine = np.eye(4) * 2
affine[:, -1] = [-90, -90, -72, 1]

gfield = gaussian_random_field(91, 109, 91, noise=noise, alpha=alpha,
                               normalize=normalize, seed=seed)
fn = make_tmpname(suffix='.nii.gz')
nib.save(nib.nifti1.Nifti1Image(gfield, affine), fn)