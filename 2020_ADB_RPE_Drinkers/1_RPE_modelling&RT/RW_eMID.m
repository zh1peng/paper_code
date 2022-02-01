function [all_RWinfo, left_ev, left_pe]=RW_eMID(EEG, retriggerEEG, eta)
% 
% caculate the RPE info - EV/PE for the EEG structure.
% RW reinforcement model to estimate trial-by-trial subjective expected
% value (EV) and prediction error (PE)
% As preprocessing and qc remove trials, EV and PE are estiamted for all
% trials and intersected with left trials. The uridx is the key used for
% this purpose.
% Input:
%   EEG: eMID EEG structure [ this function only suits the lab's data]
%Optional:
%   retriggerEEG: if urcodes are matched incorrectly.
%   eta: learning rate (default 0.3)
%
% How to use:
% EEG=pop_loadset('W:\64_EEG\EEG_data\Preprocessed\EMID\UCD_Project\DA130116BD\qc_Finalretrig_filter_DA130116BD.set')
% retriggerEEG=pop_loadset('W:\64_EEG\EEG_data\Preprocessed\EMID\UCD_Project\DA130116BD\retrig_filter_DA130116BD.set')
% [all_RWinfo, left_ev, left_pe]=RW_eMID(EEG,retriggerEEG,eta)
% author: Zhipeng
% 2018/07/18 
% 2018/07/19 user defined EEG set containing urevent is supported.
% As the retrigger removed some events but keeps the urevent in the column,
% which makes the EEG.event are not matched with EEG.urevent by the urevent
% number. Here it uses retrigger file's event as urevent.
if ~exist('retriggerEEG', 'var')
[urtype,~]=struct2vector(EEG.urevent);
else
[urtype,~]=struct2vector(retriggerEEG.event);
end

if ~exist('eta','var')
    eta=0.3;
end
cue_uridx=find(urtype==101|urtype==102|urtype==103);
fb_uridx=find(urtype==13 |urtype==23 |urtype==33 |urtype==16 |urtype==26|urtype== 36);

if length(cue_uridx)~=144
    fb_uridx=fb_uridx(fb_uridx>cue_uridx(1)); %remove fb_idx before the first cue.
end

remove_ur_fb_uridx=[]

if length(fb_uridx)>144 % remove fb_uridx that there is no cue between them.
    
    filtered_fb_uridx=[]
    fb_uridx=[1;fb_uridx]; %pad1
    for test_i=1:length(fb_uridx)-1
      cue_between=  intersect([fb_uridx(test_i):fb_uridx(test_i+1)],cue_uridx);
        if ~isempty(cue_between)
        filtered_fb_uridx=[filtered_fb_uridx;fb_uridx(test_i+1)];
        end
    end
    remove_ur_fb_uridx=setdiff( fb_uridx(2:end),filtered_fb_uridx); 
    %these urcodes have been put into the epoching, need to remove from
    % left_fb_urcode
    fb_uridx=filtered_fb_uridx;
end


cues=urtype(cue_uridx);
cues(cues==101)=20;
cues(cues==102)=-20;
cues(cues==103)=0;

fbs=urtype(fb_uridx);
fbs(fbs==13)=20;
fbs(fbs==16)=0;
fbs(fbs==23)=0;
fbs(fbs==26)=-20;
fbs(fbs==33)=0;
fbs(fbs==36)=0;



pg1=0.5; %initial p_gain
for trial_i=1:length(cues)
    ev(trial_i,1)=cues(trial_i)*pg1;
    pe(trial_i,1)=fbs(trial_i)-ev(trial_i);
    if cues(trial_i)==0
        pg1=pg1; %remains the same if it is a neutral trial
    else
        pg2=max(min(pg1+eta*(pe(trial_i)./cues(trial_i)),1.0),0.0); %make sure it is between 1 and 0
        pg1=pg2;
    end
end
all_T=[urtype(cue_uridx),cue_uridx,  urtype(fb_uridx),fb_uridx,cues, fbs,ev, pe];
all_RWinfo=array2table(all_T,'VariableNames',{'cue_mark','cue_uridx', 'fb_mark', 'fb_uridx','cue_value','fb_value','ev','pe'});
%% Get unique epoch_num

[type,latency,urcode,epoch_num]=struct2vector(EEG.event);
[~, cue_left_idx]=pop_selectevent( EEG, 'type', [101,102,103],'latency','-50<=50');
cue_left_urcode=urcode(cue_left_idx);
cue_left_epoch_num=epoch_num(cue_left_idx);
cue_left_type=type(cue_left_idx);

[~,~,ev_idx]=intersect(cue_left_urcode,cue_uridx);

ev_left=ev(ev_idx);
cue_value_left=cues(ev_idx);

left_ev_T=[cue_left_type, cue_left_urcode,cue_value_left, ev_left,cue_left_epoch_num];
left_ev=array2table(left_ev_T,'VariableNames',{'left_cue_mark','left_cue_uridx', 'left_cue_value','left_ev','ev_epoch_num'});



%     pop_select( EEG,'trial',[13,23,33,16,26,36]);
%     pop_selectevent( EEG,'type',[13,23,33,16,26,36]) latency to select correct trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  why not match them using urcode?. 
%  There will be multiple same urcode for one event, but you don't know
%  which is for the epoch!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,fb_left_idx]=pop_selectevent( EEG,'type',[13,23,33,16,26,36],'latency','-50<=50');
fb_left_urcode=urcode(fb_left_idx);
if ~isempty(remove_ur_fb_uridx)
    [~,remove_idx]=intersect(fb_left_urcode,remove_ur_fb_uridx);
    fb_left_urcode(remove_idx)=[];
    fb_left_idx(remove_idx)=[];
end
fb_left_epoch_num=epoch_num(fb_left_idx);
fb_left_type=type(fb_left_idx);

[~,~,pe_idx]=intersect(fb_left_urcode,fb_uridx);
pe_left=pe(pe_idx);
fb_value_left=fbs(pe_idx);

left_pe_T=[fb_left_type, fb_left_urcode, fb_value_left,pe_left,fb_left_epoch_num];
left_pe=array2table(left_pe_T,'VariableNames',{'left_fb_mark', 'left_fb_uridx', 'left_fb_value','left_pe','pe_epoch_num'});
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