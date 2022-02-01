function [all_results,left_results]=RT_eMID(EEG,retriggerEEG)
%% caculate the RT time for eMID EEG structure
% As preprocessing and qc remove trials, RTs are caculated for all
% trials and intersected with left trials. The uridx is the key used for
% this purpose.
%Optional:
%   retriggerEEG: if urcodes are matched incorrectly.
% 2018/07/18
% 2018/07/19 user defined EEG set containing urevent is supported.
% As the retrigger removed some events but keeps the urevent in the column,
% which makes the EEG.event are not matched with EEG.urevent by the urevent
% number. Here it uses retrigger file's event as urevent.
if ~exist('retriggerEEG', 'var')
[urtype,urlatency]=struct2vector(EEG.urevent);
else
[urtype,urlatency]=struct2vector(retriggerEEG.event);
end
%% urevent caculate all RTs.
cue_uridx=find(urtype==101|urtype==102|urtype==103);
cue_type=urtype(cue_uridx);
all_T=[cue_type,cue_uridx];

for i=1:length(cue_uridx)
    cue_lat=urlatency(cue_uridx(i));
    upper_lat=urlatency(cue_uridx(i))+1500; %search within 1500 latency interval (about 3 s)
    search_idx=find(urlatency>cue_lat&urlatency<upper_lat);
    if ~exist('retriggerEEG', 'var') %based on the urevent
        results_idx=search_idx(find(urtype(search_idx)==50)); % find target idx
        if urtype(results_idx+1)~=1 && urtype(results_idx+1)~=2 %No response (1/2) following the target
            RT(i,1)=NaN;
        else
            RT(i,1)=(urlatency(results_idx+1)-urlatency(results_idx))./0.512;
        end
        
    else % based on retriggerEEG, the markers are changed!
        results_idx=search_idx(find(urtype(search_idx)==51|urtype(search_idx)==52|urtype(search_idx)==53)); % find target idx
        if isempty(results_idx) | length(results_idx)>1
%             cue_uridx(i)=[-999];
            RT(i,1)=-999;
            continue
        end
        if urtype(results_idx+1)~=61 && urtype(results_idx+1)~=62 && urtype(results_idx+1)~=63 %No response (1/2) following the target
            RT(i,1)=NaN;
%             if urtype(results_idx+2)==61 || urtype(results_idx+2)==62
%                 RT(i,1)=(urlatency(results_idx+2)-urlatency(results_idx))./0.512;
%             end
        else
            RT(i,1)=(urlatency(results_idx+1)-urlatency(results_idx))./0.512;
        end
    end
    end
all_T=[all_T,RT];
all_results=array2table(all_T, 'VariableNames',{'cue_mark','cue_uridx','RT'});




%% intersect with left trials
[type,latency,urcode,epoch_num]=struct2vector(EEG.event);
[~, cue_left_idx]=pop_selectevent( EEG, 'type', [101,102,103],'latency','-50<=50');
cue_left_urcode=urcode(cue_left_idx);
cue_left_epoch_num=epoch_num(cue_left_idx);
cue_left_type=type(cue_left_idx);
[~,RTidx]=intersect(cue_uridx, cue_left_urcode);
left_RT=RT(RTidx);
left_T=[cue_left_type, cue_left_urcode, left_RT, cue_left_epoch_num];
left_results=array2table(left_T, 'VariableNames',{'left_cue_mark','left_cue_uridx','left_RT','left_cue_epoch_num'});
end


function varargout=struct2vector(str_in)
%example: [type, lat, urevent]=struct2vector(EEG.event);
field_names=fieldnames(str_in);
for i=1:length(field_names)
for n=1:length(str_in)
    eval(sprintf('%s(n)=str_in(n).%s;' ,field_names{i},field_names{i}));
end
end
for ii=1:nargout
varargout{ii}=eval(sprintf('transpose(%s);',field_names{ii}));
end
end