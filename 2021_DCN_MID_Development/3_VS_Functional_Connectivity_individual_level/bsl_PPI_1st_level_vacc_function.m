function []=bsl_PPI_1st_level_vacc_function(subidx)
%% 0. Setup variables
addpath(genpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12'));
addpath('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_1st_level/BSL_code')
spm_jobman('initcfg');
spm('defaults','fmri');
% warning('off','MATLAB:MKDIR:DirectoryExists');
subidfile='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_1st_level/BSL_code/subid_bsl_good_fc_file_1491.mat';
load(subidfile)
data_path= '/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_1st_level';
EPI_data_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/MID_BSL_EPI'
movement_filepath='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_info';
motion_parameter=csvread(fullfile(movement_filepath,subid{subidx},['nuisance_',subid{subidx},'.csv']));
result_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_PPI_1st_level';
if ~exist(result_path)
    mkdir(result_path)
end

% My_VOI = {...
%     'VS1_L';...
%     'VS1_R';...
%     'Postcentral_L';...
%     'Postcentral_R'};
% 
% VOI_xyz  = {...
%     [-12;14;-8];...
%     [12;14;-8];...
%     [-45;-19;55];...
%     [45;-19;55]};


My_VOI = {...
    'VS1_L';...
    'VS1_R'};

VOI_xyz  = {...
    [-12;14;-8];...
    [12;14;-8]};


%     [-12;20;-2]; vs2_L
%     [12;20;-2]; vs2_R
% %check the coordinates again
%      for n=1:length(MY_xyz)
%      disp( My_VOI{n})
%      disp(MY_xyz{n})
%      end
%% 1. F contrast (ignore) as it was done already
%% 2. extract ROI time serier adjust F contrast (effect of intrests)
cd(fullfile(data_path,subid{subidx})); load('SPM.mat');
comp={'big-no'}
cond=[1 3]
weight=[-1,1]
ppiflag = 'psychophysiologic interaction';
    for Vi=1:length(My_VOI)
        %spm_regions(SPM, my_xyz, my_name, my_con_nr+1, my_radius)
        my_spm_regions(SPM, VOI_xyz {Vi}, My_VOI{Vi}, 2, 5); %F is 1+1=2
        voi = ['VOI_' My_VOI{Vi} '_1.mat'];
        name = [My_VOI{Vi} '_' comp{1}]; %Only one comp here, no loop
        PPI = my_spm_peb_ppi(SPM, ppiflag, voi, cond, weight, name);
        clear PPI
        close all
    end
clear SPM

%% 3. Set up PPI GLM
 for Vi=1:length(My_VOI)
    
    PPI_name = ['PPI_', My_VOI{Vi} '_' comp{1}]; 
    mkdir(fullfile(result_path,subid{subidx},PPI_name));
    load(fullfile(data_path,subid{subidx},strcat(PPI_name,'.mat')));
 
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(result_path,subid{subidx},PPI_name));
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.2;
    f = spm_select('FPList', fullfile(EPI_data_path,subid{subidx}), '^swau.*\.nii$');
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(f);
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = ['PPI_interation'];
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'Pysiol_Y';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Pschol_P';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;
%     matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_parameter};
    for reg_i=1:size(motion_parameter,2)
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(reg_i+3).name=strcat('nuisance_',num2str(reg_i)); % start with 4!!!!!!!!
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(reg_i+3).val = motion_parameter(:,reg_i);
    end
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    strr=[ My_VOI{Vi},' for NO.' num2str(subidx) ' is done'];
    disp(strr)
    end
    strr=['PPI GLMs for NO.' num2str(subidx) ' are done'];
    disp(strr)


%% Estimation and Contrast
    for Vi=1:length(My_VOI)
	PPI_name = ['PPI_', My_VOI{Vi} '_' comp{1}]
    matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile(result_path,subid{subidx},PPI_name,'SPM.mat');
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{2}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = 'PPI-Interaction';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{2}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);
    clear matlabbatch
    strr=['Estimation for and Contrast of ', My_VOI{Vi}, ' for NO.', num2str(subidx), ' is done'];
    disp(strr)
    end
    strr=['Estimation for and Contrast for NO.', num2str(subidx), ' are done'];
    disp(strr)

end
