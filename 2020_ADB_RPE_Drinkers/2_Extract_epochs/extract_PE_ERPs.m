% RW PE trial extraction
% load RPE model information
load('RW_eMID_info.mat')

% load subid (cell)
load('subid.mat')
% load subject filename and filepath (cell)
load('filename_info.mat')
%% get sub RW info
for subi=1:length(subid)
    subpe=eval(['table2array(all_eMID.',subid{subi},'.pe);']);
    pos_pe_info=subpe(subpe(:,4)>0,[1,3,4,5]);
    pos_pe_value=pos_pe_info(:,3);
    pos_pe_epoch_idx=pos_pe_info(:,4);
    
    neg_pe_info=subpe(subpe(:,4)<0,[1,3,4,5]);
    neg_pe_value=neg_pe_info(:,3);
    neg_pe_epoch_idx=neg_pe_info(:,4);
    
    pos_info2save=array2table(pos_pe_info,'VariableNames',{'marker','outcome','pe_value','epoch_num'});
    neg_info2save=array2table(neg_pe_info,'VariableNames',{'marker','outcome','pe_value','epoch_num'});
    
    ORIG = pop_loadset(file_name{subi,1},file_path{subi,1});
    ORIG = eeg_checkset( ORIG );
    ORIG = pop_eegfiltnew(ORIG, [],16); % 30Hz low-pass
    ORIG=eeg_checkset(ORIG)
    ORIG= pop_reref(ORIG, [69 70] );
    ORIG = eeg_checkset( ORIG );
    
    pos_EEG= pop_select( ORIG,'time',[-0.2 2] ,'trial',pos_pe_epoch_idx);
    pos_erps=pos_EEG.data;
    pos_mean_erps=mean(pos_erps,3);
    
    neg_EEG= pop_select( ORIG,'time',[-0.2 2] ,'trial',neg_pe_epoch_idx);
    neg_erps=neg_EEG.data;
    neg_mean_erps=mean(neg_erps,3);
    
   eval('all_eMID_ERP.pos_PE(:,:,subi)=pos_mean_erps;');
   eval('all_eMID_ERP.neg_PE(:,:,subi)=neg_mean_erps;');

    eval(['all_eMID_ERP.',subid{subi},'.pos_erps=pos_erps;']);
    eval(['all_eMID_ERP.',subid{subi},'.pos_pe_value=pos_pe_value;']);
    eval(['all_eMID_ERP.',subid{subi},'.pos_pe_info=pos_info2save;']);
    
    eval(['all_eMID_ERP.',subid{subi},'.neg_erps=neg_erps;']);
    eval(['all_eMID_ERP.',subid{subi},'.neg_pe_value=neg_pe_value;']);
    eval(['all_eMID_ERP.',subid{subi},'.neg_pe_info=neg_info2save;']);

pos_EEG=[];
neg_EEG=[];
pos_erps=[];
neg_erps=[];
pos_mean_erps=[];
neg_mean_erps=[];
pos_pe_value=[];
neg_pe_value=[];
ORIG=[];
pos_info2save=[];
neg_info2save=[];

end
