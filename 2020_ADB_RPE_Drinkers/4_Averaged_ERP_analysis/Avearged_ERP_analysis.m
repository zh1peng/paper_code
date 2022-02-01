
%% load stuff needed
load('RW_eMID_PE_ERPs.mat')
pos_PE=all_eMID_ERP.pos_PE;
neg_PE=all_eMID_ERP.neg_PE;
load('chanlocs.mat')
%% load subid 
load('subid.mat')

chanlocs(65:end)=[];

color_bg=[240,249,232]./255;
cor1=[204,235,197]./255;
cor2=[168,221,181]./255;
cor3=[123,204,196]./255;
cor4=[67,162,202]./255;
cor5=[8,104,172]./255;


%% plot pos/neg PE for the whole sample
marks={'pos_PE', 'neg_PE'}
 for marki=1:length(marks)
    eval(sprintf('mean_%s=squeeze(mean(%s(:,:,:),3));',marks{marki},marks{marki}));
 end
 % select window
mean_pos_PE=mean_pos_PE(:,1:615);
mean_neg_PE=mean_neg_PE(:,1:615);
x=linspace(-200,1000,615) % x is in ms space
chans=[12 48 49 19 32 56]
chans_name=' ';
for chan=chans
    if chan==chans(1)
        chans_name=strcat(chans_name,chanlocs(chan).labels);
    end
chans_name=strcat(chans_name,',',chanlocs(chan).labels);
end
% x=1:1125;
x1=mean(mean_pos_PE(chans,:),1);
x2=mean(mean_neg_PE(chans,:),1);

w2plot1=[230:259];
w2plot2=[385:600];
figure;
% set(gcf,'color',color_bg)
hold on
column1=area(w2plot1,6*ones(size(w2plot1)),'LineStyle','none')
column1(1).FaceColor=cor1;

column2=area(w2plot2,6*ones(size(w2plot2)),'LineStyle','none')
column2(1).FaceColor=cor1;

h_line1=plot(x, x2,'color',cor4, 'LineWidth', 3)
h_line2=plot(x, x1,'color',cor2, 'LineWidth', 3)
ylim([-4,15])
% set(gca,'xtick',[-200:100:2000],'xticklabel',[-200:100:2000])
% set(gca,'XAxisLocation','origin')
% set(gca,'YAxisLocation','origin')
set(gca,'YDir','reverse')   
set(gca,'FontSize', 14, 'FontName','Arial')
% set(gca,'color',color_bg)
h_legend=legend([h_line1,h_line2],'Negative RPE','Positive RPE')
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.75,0.82,0.14,0.1])
% title([' Positive vs. Negative RPE (',chans_name,')'])
title(['Positive vs. Negative RPE'],'FontSize', 14, 'FontName','Arial')
grid on
grid minor
fig = gcf;
fig.InvertHardcopy = 'off';
saveas(gcf,'fig1.png')

%% plot binge vs. non-binge seperately for pos/neg PE
% binge_idx is the index for the binge subjects in the list
% control_idx is the index for the control subjects in the list

x=linspace(-200,1000,615);
binge_pos_PE=pos_PE(:,1:615,binge_idx);
control_pos_PE=pos_PE(:,1:615,control_idx);
binge_neg_PE=neg_PE(:,1:615,binge_idx);
control_neg_PE=neg_PE(:,1:615,control_idx);   
marks={'binge_pos_PE', 'binge_neg_PE','control_pos_PE', 'control_neg_PE'};
 for marki=1:length(marks)
    eval(sprintf('mean_%s=squeeze(mean(%s,3));',marks{marki},marks{marki}));
 end
x1=mean(mean_binge_pos_PE(chans,:),1);
x2=mean(mean_binge_neg_PE(chans,:),1);
x3=mean(mean_control_pos_PE(chans,:),1);
x4=mean(mean_control_neg_PE(chans,:),1);

figure;
set(gcf,'color',color_bg)
hold on
column1=area(w2plot1,6*ones(size(w2plot1)),'LineStyle','none')
column1(1).FaceColor=cor1;

column2=area(w2plot2,6*ones(size(w2plot2)),'LineStyle','none')
column2(1).FaceColor=cor1;


h_line2=plot(x, x2,'color',cor4, 'LineWidth', 3)
h_line4=plot(x, x4,'color',cor2, 'LineWidth', 3)
ylim([-4,15])

% set(gca,'xtick',[-200:100:2000],'xticklabel',[-200:100:2000])
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
set(gca,'YDir','reverse')   
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',color_bg)
h_legend=legend([h_line2, h_line4],'HA','LA')
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.75,0.82,0.14,0.1])
% title(['Negative RPE: BD vs. Non-BD (',chans_name,')'])
title(['Negative RPE: HA vs. LA'])
grid on
grid minor
% fig = gcf;
% fig.InvertHardcopy = 'off';
% saveas(gcf,'fig2.png')



figure;
set(gcf,'color',color_bg)
hold on
column1=area(w2plot1,6*ones(size(w2plot1)),'LineStyle','none')
column1(1).FaceColor=cor1;
column2=area(w2plot2,6*ones(size(w2plot2)),'LineStyle','none')
column2(1).FaceColor=cor1;

h_line1=plot(x, x1,'color',cor4, 'LineWidth', 3)
h_line3=plot(x, x3,'color',cor2, 'LineWidth', 3)
ylim([-4,15])

% set(gca,'xtick',[-200:100:2000],'xticklabel',[-200:100:2000])
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
set(gca,'YDir','reverse')   
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',color_bg)
h_legend=legend([h_line1, h_line3],'HA','LA')
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.75,0.82,0.14,0.1])
% title(['Positive RPE: BD vs. Non-BD (',chans_name,')'])
title(['Positive RPE: HA vs. LA'])
grid on
grid minor
% fig = gcf;
% fig.InvertHardcopy = 'off';
% saveas(gcf,'fig3.png')

%% Topograph
%% Top graph for all
diff_PE=mean_pos_PE-mean_neg_PE;
figure;
set(gcf,'color',[240,249,232]./255)
subplot(2,3,1)
topoplot(squeeze(mean(mean_pos_PE(:,220:235),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Postive RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])
subplot(2,3,2)
topoplot(squeeze(mean(mean_neg_PE(:,220:235),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Negative RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,3)
topoplot(squeeze(mean(diff_PE(:,220:235),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Postive-Negative RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])

subplot(2,3,4)
topoplot(squeeze(mean(mean_pos_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Postive RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,5)
topoplot(squeeze(mean(mean_neg_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Negative RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,6)
topoplot(squeeze(mean(diff_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Postive-Negative RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])


%% Topograph for BD and Non-BD
figure;
subplot(2,3,1)
topoplot(squeeze(mean(mean_control_pos_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD Positive RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar();
caxis([-5,5])

subplot(2,3,2)
topoplot(squeeze(mean(mean_binge_pos_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('BD Positive RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,3)
diff_PE=mean_control_pos_PE-mean_binge_pos_PE
topoplot(squeeze(mean(diff_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD minus BD Positive RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])

subplot(2,3,4)
topoplot(squeeze(mean(mean_control_pos_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD Positive RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,5)
topoplot(squeeze(mean(mean_binge_pos_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('BD Positive RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])
subplot(2,3,6)
diff_PE=mean_control_pos_PE-mean_binge_pos_PE
topoplot(squeeze(mean(diff_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD minus BD Positive RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])



figure;
set(gcf,'color',color_bg)

subplot(2,3,1)
topoplot(squeeze(mean(mean_control_neg_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD Negative RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,2)
topoplot(squeeze(mean(mean_binge_neg_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('BD Negative RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,3)
diff_PE=mean_control_neg_PE-mean_binge_neg_PE
topoplot(squeeze(mean(diff_PE(:,220:234),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD minus BD Negative RPE (230-259 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])


subplot(2,3,4)
topoplot(squeeze(mean(mean_control_neg_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD Negative RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,5)
topoplot(squeeze(mean(mean_binge_neg_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('BD Negative RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
% colorbar()
caxis([-5,5])

subplot(2,3,6)
diff_PE=mean_control_neg_PE-mean_binge_neg_PE
topoplot(squeeze(mean(diff_PE(:,299:410),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','k'});
title('Non-BD minus BD Negative RPE (384-600 ms)','FontSize', 14, 'FontName','Arial')
colorbar()
caxis([-5,5])

%% permutation results
for time_i=1:size(beta1,1)
    [~,p,~,stats1]=ttest(beta1(time_i,:));
    t1(time_i)=stats1.tstat;
    p1(time_i)=p;
    [~,~,~,stats2]=ttest(beta2(time_i,:));
    t2(time_i)=stats2.tstat;
    p2(time_i)=p;
    [~,~,~,stats3]=ttest(beta2_2(time_i,:)); 
    t2_2(time_i)=stats3.tstat;
    p2_2(time_i)=p;
end

% FDR p value
figure;
set(gcf,'color',color_bg)
hold on
x=linspace(1,600,308)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p1,0.05);
p5_line=0.05*ones(1,308);
plot(x, adj_p,'color',cor4, 'LineWidth', 2)
plot(x, p5_line,'color',cor2, 'LineWidth', 3)
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',color_bg)
h_legend=legend('FDR corrected p value','p=0.05')
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.75,0.82,0.14,0.1])
title('RPE modulation effect','FontSize', 14, 'FontName','Arial')
grid on
grid minor


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
%% additional stats on consecutive point 
a=null_t1>1.96 | null_t1 <-1.96
for per_i=1:1000
    f = find(diff([0,a(:,per_i)',0]==1));
    p = f(1:2:end-1);  % Start indices
    y = f(2:2:end)-p;  % Consecutive zeros¡¯ counts
    disp(num2str(y))
end
%% Figure 1
% Figure 1-1
x=1:308;
t196_neg=-1.96*ones(1,308);
t196_pos=1.96*ones(1,308);
figure;
set(gcf,'color',[240,249,232]./255)
hold on

ci=patch([x fliplr(x)], [up_t1 fliplr(bot_t1)],cor1)
set(ci, 'FaceColor',cor1,'EdgeColor', cor1)
% h_line1=plot(x, up_t1,'color',cor1, 'LineWidth', 2)
% h_line2=plot(x, bot_t1, 'color',cor1,'LineWidth', 2)
h_line3=plot(x, t1,'color',cor4,'LineWidth', 3);
plot(x, t196_pos,'k-.','LineWidth', 1);
plot(x, t196_neg,'k-.','LineWidth', 1);

ylim([-8,5])
title('RPE modulation effect','FontSize', 14, 'FontName','Arial')
ylabel(' {\itt} value','FontSize', 14, 'FontName','Arial')
xlabel('time (ms)','FontSize', 14, 'FontName','Arial')
h_legend=legend('95% Null t values', 'True t values', '95% Confidence Interval')
set(gca,'xtick',[0:51.2:length(x)],'xticklabel',[0:100:600])
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',[240,249,232]./255)
set(gca,'color',color_bg)
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.66,0.77,0.14,0.1])
grid on
grid minor


%Figure 1-2 
x=1:308;
figure;
set(gcf,'color',[240,249,232]./255)
hold on

ci=patch([x fliplr(x)], [up_t2 fliplr(bot_t2)],cor1)
set(ci, 'FaceColor',cor1,'EdgeColor', cor1)
% h_line1=plot(x, up_t1,'color',cor1, 'LineWidth', 2)
% h_line2=plot(x, bot_t1, 'color',cor1,'LineWidth', 2)
h_line3=plot(x, t2,'color',cor4,'LineWidth', 3);

ylim([-4,6])
title('|RPE| modulation effect','FontSize', 14, 'FontName','Arial')
ylabel(' {\itt} value','FontSize', 14, 'FontName','Arial')
xlabel('time (ms)','FontSize', 14, 'FontName','Arial')
h_legend=legend('95% Confidence Interval', 'True t values')
set(gca,'xtick',[0:51.2:length(x)],'xticklabel',[0:100:600])
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',[240,249,232]./255)
set(gca,'color',color_bg)
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.66,0.77,0.14,0.1])
grid on
grid minor
%Figure 1-3
x=1:308;
figure;
set(gcf,'color',[240,249,232]./255)
hold on

ci=patch([x fliplr(x)], [up_t2_2 fliplr(bot_t2_2)],cor1)
set(ci, 'FaceColor',cor1,'EdgeColor', cor1)
% h_line1=plot(x, up_t1,'color',cor1, 'LineWidth', 2)
% h_line2=plot(x, bot_t1, 'color',cor1,'LineWidth', 2)
h_line3=plot(x, t2_2,'color',cor4,'LineWidth', 3);

ylim([-8,4])
title('sign of RPE modulation effect','FontSize', 14, 'FontName','Arial')
ylabel(' {\itt} value','FontSize', 14, 'FontName','Arial')
xlabel('time (ms)','FontSize', 14, 'FontName','Arial')
h_legend=legend('95% Confidence Interval', 'True t values')
set(gca,'xtick',[0:51.2:length(x)],'xticklabel',[0:100:600])
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',[240,249,232]./255)
set(gca,'color',color_bg)
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.66,0.77,0.14,0.1])
grid on
grid minor


%% ERP results

tmp=1:length(t1)
frn_points=tmp(t1>up_t1);
disp(['conscutive points number: ' num2str(length(frn_points(1:end)))])
disp(['time window left: ' num2str(frn_points(1)./0.512), ' ms in idx space: ', num2str(frn_points(1)+102)])
disp(['time window right: ' num2str(frn_points(end)./0.512),' ms in idx space: ', num2str(frn_points(end)+102)])

p3_points=tmp(t1<bot_t1);
disp(['conscutive points number: ' num2str(length(p3_points(2:end)))])
disp(['time window left: ' num2str(p3_points(2)./0.512), ' ms in idx space: ', num2str(p3_points(2)+102)])
disp(['time window right: ' num2str(p3_points(end)./0.512),' ms in idx space: ', num2str(p3_points(end)+102)])


chans=[12  48 49 19 32 56] %47
g1_idx=binge_idx;
g2_idx=control_idx;

% export to SPSS
method='mean'
time_window=[220:234]; %230--259 ms 15
name_prefix='frn_pos';
data=pos_PE;
[~,pos_frn_T, ~]= erp2spss(data, chans,chanlocs,time_window, name_prefix, method,g1_idx,g2_idx);
name_prefix='frn_neg';
data=neg_PE;
[~,neg_frn_T, ~]= erp2spss(data, chans, chanlocs,time_window, name_prefix, method,g1_idx,g2_idx);

method='mean'
time_window=[299:410]; %384-600
name_prefix='p3_pos';
data=pos_PE;
[~,pos_p3_T, ~]= erp2spss(data, chans,chanlocs,time_window, name_prefix, method,g1_idx,g2_idx);
name_prefix='p3_neg';
data=neg_PE;
[~,neg_p3_T, group_T]= erp2spss(data, chans, chanlocs,time_window, name_prefix, method,g1_idx,g2_idx);
all_T=[pos_frn_T, neg_frn_T, pos_p3_T, neg_p3_T,group_T];
all_T.gender=cell2mat(cov(:,4));
writetable(all_T,'binge22_vs_control22_lp30_noref_mean_perm.csv')

