% search_dir='W:\IMAGEN_MID_data\BL_info'
% [path,filename]=filesearch_regexp(search_dir,'^mid\w*.csv');
clear;clc
% load('subid_963.mat')
load('subid_fu2_good_fc_file_1365.mat')
% subid=pad0(subid);
search_dir='H:\NAS_MRI\FU2_info';
path=fullfile(search_dir,subid);
filename=strcat('mid_',subid,'.csv');

%% import.
for n = 1: length(filename)
    file2import = fullfile(path{n},filename{n});
    delimiter = '\t'; startRow = 3;
    try
        formatSpec = '%d%s%f%f%s%f%f%f%f%s%f%f%s%f%f%f%f%[^\n\r]';
        fileID = fopen(file2import,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);
        dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]) = cellfun(@(x) num2cell(x), dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]), 'UniformOutput', false);
        tmp = [dataArray{1:end-1}];
        clearvars delimiter startRow formatSpec fileID dataArray ans;
    catch
        try
            formatSpec = '%s%s%f%f%s%f%f%f%f%s%f%f%s%f%f%f%f%[^\n\r]';
            fileID = fopen(file2import,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
            fclose(fileID);
            dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]) = cellfun(@(x) num2cell(x), dataArray([1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17]), 'UniformOutput', false);
            tmp = [dataArray{1:end-1}];
            clearvars delimiter startRow formatSpec fileID dataArray ans;
            for tmp_n=1:size(tmp,1)
                try
                    test=tmp{tmp_n,1}{1};
                    tmp{tmp_n,1}=str2double(test(regexp(test,'[\d]')));
                catch
                    tmp{tmp_n,1}=NaN;
                    
                end
            end
        end
    end
    
    
    disp(sprintf('extracting....%d / %d',n,length(filename)));
    r=size(tmp,1);
    c=size(tmp,2);
    max_trial_n=max(cell2mat(tmp(:,1)));
    csvsize(n,1)=r+c; %stardard one should be 42x17
    if  max_trial_n==42&&csvsize(n,1)>=59&c==17
        if r==43 %if 43 rows, remove last one
            tmp(43,:)=[];
        elseif r>43
            trail_num=tmp(:,1);
            [~,cut]=unique(cell2mat(trail_num),'legacy'); % remove repeat trails [cut: last ocurrance]
            if cell2mat(trail_num(cut(1)))==0
            tmp(1:cut(2)-1,:)=[]; % remove previous trials
            elseif cell2mat(trail_num(cut(1)))==1
                tmp(1:cut(1)-1,:)=[];
            end
            if isnan(cell2mat(tmp(end,1)))|cell2mat(tmp(end,1))==0
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
% 
% exclude_idx=find(~cellfun('isempty',error_sub)==1)
% subid(exclude_idx)=[];




%% 
%test too early
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'TOO_EARLY'),tmp,'Unif',0);
fu2_too_early(n)=sum(cell2mat(a));
end


%test no response
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'NO RESPONSE'),tmp,'Unif',0);
fu2_no_response(n)=sum(cell2mat(a));
end



%test too late
for n=1:size(responseType,1)
    tmp=responseType(n,:);
a=cellfun(@(x) strfind(x,'TOO_LATE'),tmp,'Unif',0);
fu2_too_late(n)=sum(cell2mat(a));
end


%test Left and Right response

for n=1:size(responseType,1)
    tmp=responseType(n,:);
right_resp=cellfun(@(x) strfind(x,'ight'),tmp,'Unif',0);
left_resp=cellfun(@(x) strfind(x,'eft'),tmp,'Unif',0);
right_num(n)=sum(cell2mat(right_resp));
left_num(n)=sum(cell2mat(left_resp));
fu2_left_right_ratio(n)=left_num(n)/right_num(n);
end
figure('units','inch','position',[0,0,18,20]);
subplot(2,2,1)
histogram(fu2_too_early,'BinWidth',1);title('T2:Too early trial number')
xline(13.5,'--r',{'Cutoff','42/3=14'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,2)
histogram(fu2_no_response,'BinWidth',1);title('T2:No response')
xline(13.5,'--r',{'Cutoff','42/3=14'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,3)
histogram(fu2_too_late,'BinWidth',1);title('T2:Too late response')
xline(20.5,'--r',{'Cutoff','42/2=21'},'LineWidth',1) % 22 is the cutoff
subplot(2,2,4)
histogram(fu2_left_right_ratio,'BinWidth',0.1);title('T2:Left/right ratio')
saveas(gcf,'FU2_trial_hist_before_QC.tif')



%% calculate
for subi=1:size(cues)
    sub_cue=cues(subi,:);
    sub_rt=rt(subi,:);
    sub_rt(sub_rt<0)=NaN % ignore negative RT. RT was recorded once.
    fu2_large_mean(subi)= nanmean(sub_rt(sub_cue==10));
    fu2_small_mean(subi)= nanmean(sub_rt(sub_cue==2));
    fu2_no_mean(subi)= nanmean(sub_rt(sub_cue==0));
    
    % find number of condition
    fu2_large_n(subi)=sum(~isnan(sub_rt(sub_cue==10)));
    fu2_small_n(subi)=sum(~isnan(sub_rt(sub_cue==2)));
    fu2_no_n(subi)=sum(~isnan(sub_rt(sub_cue==0)));
    
    % calculate
    fu2_large_std(subi)= nanstd(sub_rt(sub_cue==10));
    fu2_small_std(subi)= nanstd(sub_rt(sub_cue==2));
    fu2_no_std(subi)= nanstd(sub_rt(sub_cue==0));
end

figure('units','inch','position',[0,0,18,20]);
subplot(2,3,1)
histogram(fu2_large_n,'BinWidth',1);title('T2:Large win (trial number)')
xline(6.5,'--r',{'Valid','14/2=7'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,2)
histogram(fu2_small_n,'BinWidth',1);title('T2:Small win (trial number) ')
xline(6.5,'--r',{'Valid','14/2=7'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,3)
histogram(fu2_no_n,'BinWidth',1);title('T2:No win (trial number)')
xline(6.5,'--r',{'Valid','14/2=7'},'LineWidth',1) % 11 is the cutoff

subplot(2,3,4)
histogram(fu2_large_mean,'BinWidth',0.05);title('T2:Large win (RT)')
subplot(2,3,5)
histogram(fu2_small_mean,'BinWidth',0.05);title('T2:Small win (RT) ')
subplot(2,3,6)
histogram(fu2_no_mean,'BinWidth',0.05);title('T2:No win (RT)')
saveas(gcf,'FU2_RT_hist_before_QC.tif')

clearvars -except fu2_large_mean fu2_small_mean fu2_no_mean fu2_large_std fu2_small_std fu2_no_std ...
    fu2_too_early fu2_no_response fu2_too_late fu2_left_right_ratio fu2_large_n fu2_small_n fu2_no_n

save('fu2_RT_summary_1365before_QC.mat')