function [data_T, mean_data_T, group_T]= erp2spss(data, chans, chanlocs,time_window, name_prefix, method,g1_idx,g2_idx)
% time_window=[204:244];
% name_prefix='frn_neg_';
% data=pos_PE;
% g1_idx=binge_idx;
% g2_idx=control_idx;
% method='mean'
%  More examples:
%  [frn_neg]= erp2spss(neg_PE, [12 48 49 19 32 56], chanlocs,[220:235], 'frn_neg', 'abs_sum_area',binge_idx,control_idx);
%  [frn_pos]= erp2spss(pos_PE, [12 48 49 19 32 56], chanlocs,[220:235], 'frn_pos', 'abs_sum_area',binge_idx,control_idx);
%  [p3_neg]= erp2spss(neg_PE, [12 48 49 19 32 56], chanlocs,[299:409], 'P3_neg', 'abs_sum_area',binge_idx,control_idx);
%  [p3_pos]= erp2spss(pos_PE, [12 48 49 19 32 56], chanlocs,[299:409], 'P3_pos', 'abs_sum_area',binge_idx,control_idx);


for chani=1:length(chans)
    chan2use=chans(chani);
    chan_name{chani}=strcat(name_prefix,'_',chanlocs(chan2use).labels);
end
if strcmp(method,'mean')==1
    g1_data=squeeze(mean(data(chans, time_window,g1_idx),2))';
    g2_data=squeeze(mean(data(chans, time_window,g2_idx),2))';
    g1_mean=mean(g1_data,2);
    g2_mean=mean(g2_data,2);
    g1_g2mean=[g1_mean; g2_mean];
    mean_data_T=array2table(g1_g2mean, 'VariableNames',{[name_prefix,'_','mean_across_ele']});
    g1_g2data=[g1_data;g2_data];
    data_T=array2table(g1_g2data,'VariableNames',chan_name);
    
    
elseif strcmp(method, 'max')==1
    g1_data=squeeze(max(data(chans, time_window,g1_idx),[],2))';
    g2_data=squeeze(max(data(chans, time_window,g2_idx),[],2))';
    g1_g2data=[g1_data;g2_data];
    data_T=array2table(g1_g2data,'VariableNames',chan_name);
    
    
elseif strcmp(method, 'min')==1
    g1_data=squeeze(min(data(chans, time_window,g1_idx),[],2))';
    g2_data=squeeze(min(data(chans, time_window,g2_idx),[],2))';
    g1_g2data=[g1_data;g2_data];
    data_T=array2table(g1_g2data,'VariableNames',chan_name);
    
    
elseif strcmp(method, 'abs_sum_area')==1
    g1_tmp_data=squeeze(mean(data(chans,:,g1_idx),1));
    g2_tmp_data=squeeze(mean(data(chans, :,g2_idx),1));
    g1_data=trapz(time_window,abs(g1_tmp_data(time_window,:)));
    g2_data=trapz(time_window,abs(g2_tmp_data(time_window,:)));
    
    g1_g2data=[g1_data,g2_data]';
    mean_data_T=array2table(g1_g2data, 'VariableNames',{[name_prefix,'_','abs_sum_area']});
    data_T=mean_data_T;
    
elseif strcmp(method, 'sum_area')==1
    g1_tmp_data=squeeze(mean(data(chans,:,g1_idx),1));
    g2_tmp_data=squeeze(mean(data(chans, :,g2_idx),1));
    g1_data=trapz(time_window,g1_tmp_data(time_window,:));
    g2_data=trapz(time_window,g2_tmp_data(time_window,:));
    
    g1_g2data=[g1_data,g2_data]';
    mean_data_T=array2table(g1_g2data, 'VariableNames',{[name_prefix,'_','sum_area']});
    data_T=mean_data_T;
end
group_T=array2table([ones(length(g1_idx),1);2*ones(length(g2_idx),1)],'VariableNames',{'group'});
end



