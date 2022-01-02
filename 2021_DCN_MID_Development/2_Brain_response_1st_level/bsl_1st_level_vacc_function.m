function []=bsl_1st_level_vacc_function(indx)
addpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12');
spm_jobman('initcfg');
spm('defaults','fmri');
% warning('off','MATLAB:MKDIR:DirectoryExists');

data_path= '/gpfs2/scratch/zcao4/IMAGEN_MID_data/MID_BSL_EPI';
movement_filepath='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_info';
result_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/BSL_1st_level';


subidfile='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/1st_level/BSL_code/subid2use_1509.mat';
behavfile='/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/1st_level/BSL_code/behav_info.mat';
load(subidfile)
load(behavfile)


if ~exist(result_path)
    mkdir(result_path)
end


%Regressors
%Subject's type and onsets of cues and feedbacks
all_subjects=allsubjects.subjects';
cue_type=allsubjects.cues;
outcome_type=allsubjects.rewards;
anticipation_onsets=allsubjects.anticipationOnsets;
outcome_onsets=allsubjects.feedbackOnsets;


substr=['mid_', subid{indx}, '.csv'];
sub_index=find(strcmp(all_subjects, substr)==1);%Find data indx corresponding to the subject

sub_cue=cue_type(sub_index,:);
sub_anticipation_onsets=anticipation_onsets(sub_index,:);
sub_outcome=outcome_type(sub_index,:);
sub_outcome_onsets=outcome_onsets(sub_index,:);

big_win_index=find(sub_cue==10);
small_win_index=find(sub_cue==2);
no_win_index=find(sub_cue==0);

%Anticipation Regressors big-small-no win cues
big_win_cue_onsets=sub_anticipation_onsets(big_win_index);
small_win_cue_onsets=sub_anticipation_onsets(small_win_index);
no_win_cue_onsets=sub_anticipation_onsets(no_win_index);

%Outcome Regressors 2x2+1(No Win)
big_win_outcome=sub_outcome(big_win_index);
big_win_outcome_onsets=sub_outcome_onsets(big_win_index);
Hit_big_win_onsets=big_win_outcome_onsets(find(big_win_outcome==10));
Mis_big_win_onsets=big_win_outcome_onsets(find(big_win_outcome==0));

small_win_outcome=sub_outcome(small_win_index);
small_win_outcome_onsets=sub_outcome_onsets(small_win_index);
Hit_small_win_onsets=small_win_outcome_onsets(find(small_win_outcome==2));
Mis_small_win_onsets=small_win_outcome_onsets(find(small_win_outcome==0));

no_win_onsets=sub_outcome_onsets(no_win_index);


% f = spm_select('FPList', fullfile(data_path,subid{indx}), '^swea.*\.nii$');
f = spm_select('FPList', fullfile(data_path,subid{indx}), '^swau.*\.nii$');
motion_parameter=csvread(fullfile(movement_filepath,subid{indx},['nuisance_',subid{indx},'.csv']));

sub_result_path=fullfile(result_path,subid{indx});

if ~exist(sub_result_path)
    mkdir(sub_result_path)
end


matlabbatch{1}.spm.stats.fmri_spec.dir =cellstr(sub_result_path);
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(f);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;


%Anticipation condition ->Durations:  4 for anticipation and 1.5 for feedback
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'no_win';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = no_win_cue_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 4;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'small_win';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = small_win_cue_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 4;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name = 'big_win';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset = big_win_cue_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration = 4;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;

%Outcome condition
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name = 'big_win_hit';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset = Hit_big_win_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration = 1.45;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name = 'big_win_mis';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset = Mis_big_win_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration = 1.45;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name = 'small_win_hit';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset = Hit_small_win_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration = 1.45;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name = 'small_win_mis';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset = Mis_small_win_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration = 1.45;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name = 'no_win_outcome';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset = no_win_onsets';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration = 1.45;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod = 0;
% matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {motion_parameter};
for reg_i=1:size(motion_parameter,2)
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(reg_i).name=strcat('nuisance_',num2str(reg_i));    
matlabbatch{1}.spm.stats.fmri_spec.sess.regress(reg_i).val = motion_parameter(:,reg_i);
end
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Effect of interst';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(8);
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'anticipation_big-no';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
try
spm_jobman('run',matlabbatch);
clearvars -except indx data_path movement_filepath result_path subid all_subjects cue_type outcome_type anticipation_onsets allsubjects.anticipationOnsets outcome_onsets
strr=['No.' num2str(indx) ' subject is done'];
disp(strr)
catch
    error_sub=subid{indx};
    save(['/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/1st_level/BSL_code/error_',num2str(indx),'.mat'], error_sub)
end
end

