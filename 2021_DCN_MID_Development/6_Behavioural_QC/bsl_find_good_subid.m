% Subject inclusion criteria:
%  BSL: Too early < 66/3=22
%          No response < 66/3=22
%          Too late < 66/2=33
%  FU2: Too early < 42/3 =14
%           No response < 42/3=14
%           Too late < 42/2=21
% BSL_condition_n>10
% FU2_condition_n>6
 load('subid_bsl_good_fc_file_1491.mat')
load('bsl_RT_summary_1491before_QC.mat')
good_response=bsl_too_early<22&bsl_no_response<22&bsl_too_late<33
good_condition=bsl_large_n>10&bsl_small_n>10&bsl_no_n>10&bsl_large_mean>0.01
bsl_good=good_response&good_condition
subid=subid(bsl_good)
save('bsl_good_behav_subid1304.mat','subid')




