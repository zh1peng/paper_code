% load epoched data based on PE values
load('RPE_epoched_data.mat')
% load channel information
load('chanlocs.mat')
chanlocs(69:end)=[];

clearvars -except all_eMID_ERP

time_window=[102:409] % 200-800 (200 baseline)
chans=[48]
for subi=1:length(subid)
    sub_pos_pe=eval(['all_eMID_ERP.',subid{subi},'.pos_pe_value;']);
    sub_neg_pe=eval(['all_eMID_ERP.',subid{subi},'.neg_pe_value;']);
    sub_pos_erps=eval(['all_eMID_ERP.',subid{subi},'.pos_erps;']);
    sub_neg_erps=eval(['all_eMID_ERP.',subid{subi},'.neg_erps;']);
    chan_sub_pos_erps=squeeze(mean(sub_pos_erps(chans,time_window,:),1));
    chan_sub_neg_erps=squeeze(mean(sub_neg_erps(chans,time_window,:),1));
    chan_erps=cat(2, chan_sub_pos_erps, chan_sub_neg_erps);
    all_pe=[sub_pos_pe; sub_neg_pe];
    
    for time_i=1:size(chan_erps,1)
        model1=[ones(length(all_pe),1),all_pe];
        model2=[ones(length(all_pe),1),abs(all_pe),sign(all_pe)];    
        tmp_beta1=pinv(model1)*chan_erps(time_i,:)';
        tmp_beta2=pinv(model2)*chan_erps(time_i,:)';
        beta1(time_i, subi)=tmp_beta1(2);
        beta2(time_i, subi)=tmp_beta2(2);
        beta2_2(time_i, subi)=tmp_beta2(3);
        for perm_i=1:1000
            shuffle_pe=all_pe(randperm(length(all_pe)));
            null_model1=[ones(length(shuffle_pe),1),shuffle_pe];
            null_model2=[ones(length(shuffle_pe),1),abs(shuffle_pe),sign(shuffle_pe)];
            null_tmp_beta1=pinv(null_model1)*chan_erps(time_i,:)';
            null_tmp_beta2=pinv( null_model2)*chan_erps(time_i,:)';
            null_beta1(time_i, perm_i, subi)=null_tmp_beta1(2);
            null_beta2(time_i, perm_i,subi)=null_tmp_beta2(2);
            null_beta2_2(time_i, perm_i,subi)=null_tmp_beta2(3);
        end
        disp(sprintf('subject: %d, at timepoint: %d', subi, time_i))
    end
end



%% Single trial analysis for whole group
perm_num=1000;
for time_i=1:size(null_beta1,1)
    for perm_i=1:perm_num
    [~,~,~,stats1]=ttest(squeeze(null_beta1(time_i,perm_i,:)));
    null_t1(time_i,perm_i)=stats1.tstat;
    [~,~,~,stats2]=ttest(squeeze(null_beta2(time_i,perm_i,:)));
    null_t2(time_i,perm_i)=stats2.tstat;
    [~,~,~,stats3]=ttest(squeeze(null_beta2_2(time_i,perm_i,:))); 
    null_t2_2(time_i,perm_i)=stats3.tstat;
    end
    up_t1(time_i)=prctile(null_t1(time_i,:), 97.5);
    bot_t1(time_i)=prctile(null_t1(time_i,:), 2.5);
    up_t2(time_i)=prctile(null_t2(time_i,:), 97.5);
    bot_t2(time_i)=prctile(null_t2(time_i,:), 2.5);
    up_t2_2(time_i)=prctile(null_t2_2(time_i,:), 97.5);
    bot_t2_2(time_i)=prctile(null_t2_2(time_i,:), 2.5);    
end




%% plot
time_window=102:409
x=1:1:length(time_window)

% http://colorbrewer2.org/#type=sequential&scheme=GnBu&n=4
figure;
set(gcf,'color',[240,249,232]./255)
hold on
h_line1=plot(x, up_t1,'color',[186,228,188]./255, 'LineWidth', 3)
h_line2=plot(x, bot_t1, 'color',[186,228,188]./255,'LineWidth', 3)
h_line3=plot(x, t1,'color',[43,140,190]./255,'LineWidth', 3);
ylim([-8,5])
title('RPE modulation effect with permuation test','FontSize', 14, 'FontName','Arial')
ylabel(' {\itt} value','FontSize', 14, 'FontName','Arial')
xlabel('time (ms)','FontSize', 14, 'FontName','Arial')
h_legend=legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:600])
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',[240,249,232]./255)
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',[240,249,232]./255,...
    'EdgeColor', [240,249,232]./255,...
    'Position',[0.66,0.77,0.14,0.1])
grid on
grid minor


figure;
plot(x, up_t2,'r',...
    x, bot_t2, 'r',...
    x, t2,'b');
ylim([-5,5])
title('|RPE| modulation effect with permuation test')
ylabel('t value')
xlabel('time (ms)')
legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:600])

figure;
plot(x, up_t2_2,'r',...
    x, bot_t2_2, 'r',...
    x, t2_2,'b');
ylim([-8,5])
title('sign RPE modulation effect with permuation test')
ylabel('t value')
xlabel('time (ms)')
legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:600])

%% single trial analysis for HA vs. LA

% subidx is the index of HA subjects in the list
binge_subidx=binge_idx; 
control_subidx=control_idx;


for time_i=1:size(beta1,1)
binge_beta1=beta1(time_i, binge_subidx);
binge_beta2= beta2(time_i, binge_subidx);
binge_beta2_2=beta2_2(time_i, binge_subidx);

control_beta1=beta1(time_i, control_subidx);
control_beta2= beta2(time_i, control_subidx);
control_beta2_2=beta2_2(time_i, control_subidx);

    [~,p1(time_i)]=ttest2(binge_beta1,control_beta1);
    [~,p2(time_i)]=ttest2(binge_beta2,control_beta2);
    [~,p2_2(time_i)]=ttest2(binge_beta2_2,control_beta2_2);

end


binge2plot=beta1(:,binge_subidx);
control2plot=beta1(:,control_subidx);
figure
hold on
for ii=1:size(binge2plot,2)
x=1:1:length(time_window)
plot(x, binge2plot(:,ii),'r')
end

test=binge2plot(140,:)


for ii=1:size(control2plot,2)
plot(x, control2plot(:,ii),'g')
end



[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p1,0.05);
x=1:1:length(time_window)
figure;
plot(x, adj_p);
 ylim([0,1])
hline(0.05,'r-')
title('RPE modulation effect')
ylabel('FDR corrected p value')
xlabel('time (ms)')
legend('FDR corrected p value')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])


[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p2,0.05);
x=1:1:length(time_window)
figure;
plot(x, adj_p);
 ylim([0,1])
hline(0.05,'r-')
title('|RPE| modulation effect')
ylabel('FDR corrected p value')
xlabel('time (ms)')
legend('FDR corrected p value')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p2_2,0.05);
x=1:1:length(time_window)
figure;
plot(x, adj_p);
 ylim([0,1])
hline(0.05,'r-')
title('sign of RPE modulation effect')
ylabel('FDR corrected p value')
xlabel('time (ms)')
legend('FDR corrected p value')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])

%% permutation approach 1: use fake betas
binge_subidx=binge_idx;
control_subidx=control_idx;
for time_i=1:size(beta1,1)
        binge_beta1=beta1(time_i, binge_subidx);
        binge_beta2= beta2(time_i, binge_subidx);
        binge_beta2_2=beta2_2(time_i, binge_subidx);
        
        control_beta1=beta1(time_i, control_subidx);
        control_beta2= beta2(time_i, control_subidx);
        control_beta2_2=beta2_2(time_i, control_subidx);
        [~,~,~,stats1]=ttest2(binge_beta1,control_beta1);
        [~,~,~,stats2]=ttest2(binge_beta2,control_beta2);
        [~,~,~,stats3]=ttest2(binge_beta2_2,control_beta2_2);
        t1(time_i)=stats1.tstat;
        t2(time_i)=stats2.tstat;
        t2_2(time_i)=stats3.tstat;
end


for time_i=1:size(beta1,1)
    for perm_i=1:perm_num
        binge_beta1=squeeze(null_beta1(time_i, perm_i,binge_subidx));
        binge_beta2= squeeze(null_beta2(time_i, perm_i,binge_subidx));
        binge_beta2_2=squeeze(null_beta2_2(time_i, perm_i,binge_subidx)); 
        control_beta1=squeeze(null_beta1(time_i, perm_i,control_subidx));
        control_beta2= squeeze(null_beta2(time_i, perm_i,control_subidx));
        control_beta2_2=squeeze(null_beta2_2(time_i, perm_i,control_subidx));
        [~,~,~,stats1]=ttest2(binge_beta1,control_beta1);
        [~,~,~,stats2]=ttest2(binge_beta2,control_beta2);
        [~,~,~,stats3]=ttest2(binge_beta2_2,control_beta2_2);
        null_t1(time_i,perm_i)=stats1.tstat;
        null_t2(time_i,perm_i)=stats2.tstat;
        null_t2_2(time_i,perm_i)=stats3.tstat;
    end
    
    up_t1(time_i)=prctile(null_t1(time_i,:), 97.5);
    bot_t1(time_i)=prctile(null_t1(time_i,:), 2.5);
    
    up_t2(time_i)=prctile(null_t2(time_i,:), 97.5);
    bot_t2(time_i)=prctile(null_t2(time_i,:), 2.5);
    
    up_t2_2(time_i)=prctile(null_t2_2(time_i,:), 97.5);
    bot_t2_2(time_i)=prctile(null_t2_2(time_i,:), 2.5);
end

time_window=102:409
x=1:1:length(time_window)
figure;
plot(x, up_t1,'r',...
    x, bot_t1, 'r',...
    x, t1,'b');
ylim([-8,5])
title('RPE modulation effect with permuation test')
ylabel('t value')
xlabel('time (ms)')
legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])


figure;
plot(x, up_t2,'r',...
    x, bot_t2, 'r',...
    x, t2,'b');
ylim([-5,5])
title('|RPE| modulation effect with permuation test')
ylabel('t value')
xlabel('time (ms)')
legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])


figure;
plot(x, up_t2_2,'r',...
    x, bot_t2_2, 'r',...
    x, t2_2,'b');
ylim([-8,5])
title('sign RPE modulation effect with permuation test')
ylabel('t value')
xlabel('time (ms)')
legend('uppder 95% boundary', 'lower 5% boundary', 'true t values')
set(gca,'xtick',[0:51.2:length(time_window)],'xticklabel',[0:100:800])

%% permutation approach 2: shuffle true betas [got similar results as approach 1]
for time_i=1:size(beta1,1)
        binge_beta1=beta1(time_i, binge_idx);
        binge_beta2= beta2(time_i, binge_idx);
        binge_beta2_2=beta2_2(time_i, binge_idx);
        control_beta1=beta1(time_i, control_idx);
        control_beta2= beta2(time_i, control_idx);
        control_beta2_2=beta2_2(time_i, control_idx);
        [~,~,~,stats1]=ttest2(binge_beta1,control_beta1);
        [~,~,~,stats2]=ttest2(binge_beta2,control_beta2);
        [~,~,~,stats3]=ttest2(binge_beta2_2,control_beta2_2);
        t1(time_i)=stats1.tstat;
        t2(time_i)=stats2.tstat;
        t2_2(time_i)=stats3.tstat;
        
        for perm_i=1:perm_num
         perm_idx=randperm(sum(binge_idx|control_idx));
         null_binge_idx=perm_idx(1:sum(binge_idx));
         null_control_idx=perm_idx(sum(control_idx):end);
        null_binge_beta1=beta1(time_i, null_binge_idx);
        null_binge_beta2= beta2(time_i, null_binge_idx);
        null_binge_beta2_2=beta2_2(time_i, null_binge_idx);
        null_control_beta1=beta1(time_i, null_control_idx);
        null_control_beta2= beta2(time_i, null_control_idx);
        null_control_beta2_2=beta2_2(time_i, null_control_idx);
        [~,~,~,stats1]=ttest2(null_binge_beta1,null_control_beta1);
        [~,~,~,stats2]=ttest2(null_binge_beta2,null_control_beta2);
        [~,~,~,stats3]=ttest2(null_binge_beta2_2,null_control_beta2_2);
        null_t1(time_i,perm_i)=stats1.tstat;
        null_t2(time_i,perm_i)=stats2.tstat;
        null_t2_2(time_i,perm_i)=stats3.tstat;
        end
        
    up_t1(time_i)=prctile(null_t1(time_i,:), 97.5);
    bot_t1(time_i)=prctile(null_t1(time_i,:), 2.5);
    
    up_t2(time_i)=prctile(null_t2(time_i,:), 97.5);
    bot_t2(time_i)=prctile(null_t2(time_i,:), 2.5);
    
    up_t2_2(time_i)=prctile(null_t2_2(time_i,:), 97.5);
    bot_t2_2(time_i)=prctile(null_t2_2(time_i,:), 2.5);

end
