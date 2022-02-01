% load subid and qced file names 
load('subid.mat')
load('filename_info.mat')

% or create from subid
% data_path='W:\64_EEG\EEG_data\Preprocessed\EMID'
% file_id=cellfun(@(x) strcat('qc_Finalretrig_filter_',x),cov(:,1),'Unif',0)
% for subi=1:length(file_id)
%     [file_path(subi,1),file_name(subi,1)]=filesearch_substring(data_path,[file_id{subi},'.set']);
% end

for n=1:length(subid)
    
    EEG=pop_loadset(file_name{n,1},file_path{n,1});
    EEG=eeg_checkset(EEG);
    retriggerEEG=pop_loadset(file_name{n,1}(9:end),file_path{n,1})

    [all_comp,ev,pe]=RW_eMID(EEG,retriggerEEG)
    eval(['all_eMID.',subid{n},'.all_comp=all_comp;'])
    eval(['all_eMID.',subid{n},'.ev=ev;'])
    eval(['all_eMID.',subid{n},'.pe=pe;'])

end
for n=1:length(subid)
    
    EEG=pop_loadset(file_name{n,1},file_path{n,1});
    EEG=eeg_checkset(EEG);
    retriggerEEG=pop_loadset(file_name{n,1}(9:end),file_path{n,1})
    [all_RT, left_RT]=RT_eMID(EEG,retriggerEEG)
    eval(['all_eMID.',subid{n},'.all_RT=all_RT;'])
    eval(['all_eMID.',subid{n},'.left_RT=left_RT;']) 
    EEG=[];
    retriggerEEG=[];
end
        