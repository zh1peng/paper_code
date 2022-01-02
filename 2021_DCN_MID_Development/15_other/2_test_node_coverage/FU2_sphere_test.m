addpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12');
data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/FU2_1st_level';
load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/sphere_coverage_test/subid_FU2_1367.mat');
load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/sphere_coverage_test/node166_mnis.mat');
vols=fullfile(data_path, subid, 'mask.nii')
result_3=sphere_testing(vols, mnis,3);
result_5=sphere_testing(vols, mnis,5);
save('FU2_coverage.mat', 'result_3','result_5')