% search_dir='W:\IMAGEN_MID_data\BL_info'
% [path,filename]=filesearch_regexp(search_dir,'^mid\w*.csv');
clc;clear
% load('subid_963.mat')
 load('bsl_good_behav_subid1304.mat')
 % subid=pad0(subid);
search_dir='H:\NAS_MRI\BL_info';
path=fullfile(search_dir,subid);
filename=strcat('mid_',subid,'.csv');

%% import.
for n = 1: length(filename)
    file2import = fullfile(path{n},filename{n});
    delimiter = '\t'; startRow = 3;
    formatSpec = '%f%s%f%f%s%f%f%f%f%s%f%f%s%f%f%f%f%[^\n\r]';
    fileID = fopen(file2import,'r');
    try
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);
        dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]) = cellfun(@(x) num2cell(x), dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]), 'UniformOutput', false);
        tmp = [dataArray{1:end-1}];
        clearvars delimiter startRow formatSpec fileID dataArray ans;
    catch
        cannot_read_csv{n}=file2import;
        continue
    end
    
    disp(sprintf('extracting....%d / %d',n,length(filename)));
    r=size(tmp,1);
    c=size(tmp,2);
    csvsize(n,1)=r+c; %stardard one should be 66x17
    if  csvsize(n,1)>=83&c==17
        if r==67 %if 67 rows, remove last one
            tmp(67,:)=[];
        elseif r>67
            trail_num=tmp(:,1);
            [~,cut]=unique(cell2mat(trail_num),'legacy'); % remove repeat trails [cut: last ocurrance]
            tmp(1:cut(1)-1,:)=[]; % remove previous trials
            if isnan(cell2mat(tmp(end,1)))
                tmp(end,:)=[];
            end
                
        end
        subjects{n}=filename{n};
        tmp_cues=tmp(:,2); %Trial Category
        for m=1:length(tmp_cues)
            if strcmp(tmp_cues{m},'BIG_WIN')==1
                cues(n,m)=10;
            elseif strcmp(tmp_cues{m},'SMALL_WIN')==1
                cues(n,m)=2;
            elseif strcmp(tmp_cues{m},'NO_WIN')==1
                cues(n,m)=0;
            end
        end
        clear tmp_cues
        rewards(n,:)=transpose(cell2mat(tmp(:,14))); %amount
        anticipationOnsets(n,:)=transpose(cell2mat(tmp(:,6))/1000); %anticipaiton phase start time
        feedbackOnsets(n,:) =transpose(cell2mat(tmp(:,12))/1000); %feedback phase start time
        targetOnsets(n,:) = transpose(cell2mat(tmp(:,8))/1000); %target phase start time
        fixationOnsets(n,:)= transpose(cell2mat(tmp(:,15))/1000);
        rt(n,:)=transpose((cell2mat(tmp(:,11))-cell2mat(tmp(:,8)))/1000);
        responseType(n,:)=transpose(tmp(:,10));
        outcomeType(n,:)=transpose(tmp(:,13));
        
%         v2save(:,1)=cues(n,:)';
%         v2save(:,2)=rewards(n,:)';
%         v2save(:,3)=anticipationOnsets(n,:)';
%         v2save(:,4)=feedbackOnsets(n,:)';
%         v_savename=fullfile('W:\IMAGEN_MID_data\BL_info',subid{n},strcat('events_',num2str(subid{n}),'.csv'));
%         csvwrite(v_savename,v2save);
        
    else
        error_sub{n}=file2import;
    end
    clear tmp r c
end

% cannot_read_csv=cannot_read_csv(~cellfun('isempty',cannot_read_csv));
% % winopen(cannot_read_csv{1})
% error_sub=error_sub(~cellfun('isempty',error_sub));


%test too early
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'TOO_EARLY'),tmp,'Unif',0);
bsl_too_early(n)=sum(cell2mat(a));
end


%test no response
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'NO RESPONSE'),tmp,'Unif',0);
bsl_no_response(n)=sum(cell2mat(a));
end



%test too late
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'TOO_LATE'),tmp,'Unif',0);
bsl_too_late(n)=sum(cell2mat(a));
end


%test Left and Right response

for n=1:size(responseType,1)
    tmp=responseType(n,:);
right_resp=cellfun(@(x) strfind(x,'ight'),tmp,'Unif',0);
left_resp=cellfun(@(x) strfind(x,'eft'),tmp,'Unif',0);
right_num(n)=sum(cell2mat(right_resp));
left_num(n)=sum(cell2mat(left_resp));
bsl_left_right_ratio(n)=left_num(n)/right_num(n);
end
figure('units','inch','position',[0,0,18,20]);
subplot(2,2,1)
histogram(bsl_too_early,'BinWidth',1);title('T1:Too early trial number')
xline(21.5,'--r',{'Cutoff','66/3=22'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,2)
histogram(bsl_no_response,'BinWidth',1);title('T1:No response')
xline(21.5,'--r',{'Cutoff','66/3=22'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,3)
histogram(bsl_too_late,'BinWidth',1);title('T1:Too late response')
xline(32.5,'--r',{'Cutoff','66/2=33'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,4)
histogram(bsl_left_right_ratio,'BinWidth',0.1);title('T1:Left/right ratio')
saveas(gcf,'BSL_trial_hist_after_QC.tif')

%% calculate
for subi=1:size(cues)
    sub_cue=cues(subi,:);
    sub_rt=rt(subi,:);
    sub_rt(sub_rt<0)=NaN % ignore negative RT that is response too early
    bsl_large_mean(subi)= nanmean(sub_rt(sub_cue==10));
    bsl_small_mean(subi)= nanmean(sub_rt(sub_cue==2));
    bsl_no_mean(subi)= nanmean(sub_rt(sub_cue==0));
    
    % find number of condition
    bsl_large_n(subi)=sum(~isnan(sub_rt(sub_cue==10)));
    bsl_small_n(subi)=sum(~isnan(sub_rt(sub_cue==2)));
    bsl_no_n(subi)=sum(~isnan(sub_rt(sub_cue==0)));
    
    % calculate
    bsl_large_std(subi)= nanstd(sub_rt(sub_cue==10));
    bsl_small_std(subi)= nanstd(sub_rt(sub_cue==2));
    bsl_no_std(subi)= nanstd(sub_rt(sub_cue==0));
end

figure('units','inch','position',[0,0,18,20]);
subplot(2,3,1)
histogram(bsl_large_n,'BinWidth',1);title('T1:Large win (trial number)')
xline(10.5,'--r',{'Valid','22/2=11'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,2)
histogram(bsl_small_n,'BinWidth',1);title('T1:Small win (trial number) ')
xline(10.5,'--r',{'Valid','22/2=11'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,3)
histogram(bsl_no_n,'BinWidth',1);title('T1:No win (trial number)')
xline(10.5,'--r',{'Valid','22/2=11'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,4)
histogram(bsl_large_mean,'BinWidth',0.05);title('T1:Large win (RT)')
subplot(2,3,5)
histogram(bsl_small_mean,'BinWidth',0.05);title('T1:Small win (RT) ')
subplot(2,3,6)
histogram(bsl_no_mean,'BinWidth',0.05);title('T1:No win (RT)')
saveas(gcf,'BSL_RT_hist_after_QC.tif')



clearvars -except bsl_large_mean bsl_small_mean bsl_no_mean bsl_large_std bsl_small_std bsl_no_std ...
    bsl_too_early bsl_no_response bsl_too_late bsl_left_right_ratio bsl_large_n bsl_small_n bsl_no_n
save('bsl_RT_summary_after_QC.mat')
