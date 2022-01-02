% Subject inclusion criteria:
%  BSL: Too early < 66/3=22
%          No response < 66/3=22
%          Too late < 66/2=33
%  FU2: Too early < 42/3 =14
%           No response < 42/3=14
%           Too late < 42/2=21
% BSL_condition_n>10
% FU2_condition_n>6
load('fu2_RT_summary_1365before_QC.mat')
load('subid_fu2_good_fc_file_1365.mat')
good_response=fu2_too_early<14&fu2_no_response<14&fu2_too_late<21
good_condition=fu2_large_n>6&fu2_small_n>6&fu2_no_n>6&fu2_large_mean>0.01
fu2_good=good_response&good_condition
subid=subid(fu2_good)
save('fu2_good_behav_subid1241.mat','subid')

writematrix(fu2_good','fu2_good.csv')


