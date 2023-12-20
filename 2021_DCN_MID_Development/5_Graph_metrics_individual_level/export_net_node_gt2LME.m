load('network_results.mat')
load('F:\Google Drive\post-doc\MID_network_FU2\process_code\fc_gt_analysis\network_extraction_new_bu\Nodal analysis\node166_info.mat')
var2plot={'Cp', 'Lp', 'locE', 'gE',  'deg', 'bw', 'mod'}
var2lme=strcat('a',var2plot)

result_T=table;
for vari=1:length(var2lme)
    bl_data=eval(['bsl_net.',var2lme{vari}]);
    fu_data=eval(['fu2_net.',var2lme{vari}]);
    eval(['result_T.net_',var2lme{vari},'=[transpose(bl_data); transpose(fu_data)];'])
end
for vari=1:length(var2lme)
    if strcmp(var2lme{vari},'amod')~=1
for node_i=1:166
node_bsl_data=eval(['bsl_node.',var2lme{vari},'(node_i,:)']);
node_fu2_data=eval(['fu2_node.',var2lme{vari},'(node_i,:)']);
var_label=sprintf('result_T.node_%s_node%d',var2lme{vari},node_i)
 eval([var_label,'=[transpose(node_bsl_data); transpose(node_fu2_data)];'])
end
    end
end

writetable(result_T,'all_gt_measures.csv')