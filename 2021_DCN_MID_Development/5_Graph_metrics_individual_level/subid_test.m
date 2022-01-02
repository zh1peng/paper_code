load('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/fc_analysis/bsl_fc_code/subid_bsl_good_fc_file_1491.mat');
data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_1st_level';
fc_file=fullfile(data_path, subid, 'fc_map_166.mat');
fc_test=cellfun(@exist, fc_file, 'Unif',0);
save('bsl_good_fc_files.mat','fc_test','subid')