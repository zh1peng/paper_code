% make final plots new version
color_bg=[240,249,232]./255;
cor1=[204,235,197]./255;
cor2=[168,221,181]./255;
cor3=[123,204,196]./255;
cor4=[67,162,202]./255;
cor5=[8,104,172]./255;


color_bg=[1,1,1]
%% Fig1. ERP Grad average for timepoint selection
chans=[12 48 49 19 32 56]
chans_name=' ';
for chan=chans
    if chan==chans(1)
        chans_name=strcat(chans_name,chanlocs(chan).labels);
    end
chans_name=strcat(chans_name,',',chanlocs(chan).labels);
end
x=1:1125;
x1=mean(mean_pos_PE(chans,:),1);
x2=mean(mean_neg_PE(chans,:),1);
x_mean=(x1+x2)./2

frn_t=[204:256]
p3_t=[256:409]
[frn_value, frn_idx]=min(x_mean(:,frn_t),[],2)
frn_real_idx=frn_t(frn_idx);
[p3_value, p3_idx]=max(x_mean(:,p3_t),[],2)
p3_real_idx=p3_t(p3_idx); % idx is in point space

x=linspace(-200,2000,1125) % x is in ms space
figure;
hold on
h_line1=plot(x, x_mean,'color',cor2, 'LineWidth', 2)
plot(x(frn_real_idx),frn_value,'s','color', cor3,'MarkerSize',12)
plot(x(p3_real_idx),p3_value,'s','color', cor5,'MarkerSize',12)
ylim([-2,5])
set(gca,'xtick',[-200:100:2000],'xticklabel',[-200:100:2000])
set(gca,'XAxisLocation','origin')
set(gca,'YAxisLocation','origin')
set(gca,'YDir','reverse')   
set(gca,'FontSize', 14, 'FontName','Arial')
set(gca,'color',color_bg)
h_legend=legend('Grand Avearge','FRN Peak','P3 Peak')
set(h_legend,'FontSize', 14, 'FontName','Arial',...
    'color',color_bg,...
    'EdgeColor', color_bg,...
    'Position',[0.75,0.82,0.14,0.1])
title(['Grand average of feedback trials across ',chans_name])
grid on
text(x(frn_real_idx)+30, double(round(frn_value,1)), [num2str(round(x(frn_real_idx))),'ms'],...
                                'FontSize', 14, 'FontName','Arial','color', cor5)
text(x(p3_real_idx)+30, double(round(p3_value,1)),  [num2str(round(x(p3_real_idx))),'ms'],...
                                    'FontSize', 14, 'FontName','Arial','color', cor5)




mean2plot=(mean_pos_PE+mean_neg_PE)./2
% FRN top
figure;
subplot(1,2,1)
start=(200+200)*0.512;
endpoint=(200+300)*0.512;
topoplot(squeeze(mean(mean2plot(:,start:endpoint),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','black'});
title(sprintf('Grand average of feedback trials (%s-%s ms)', num2str(200),num2str(200+100)))
caxis([-8,8])
colorbar()
subplot(1,2,2)
start=(200+300)*0.512;
endpoint=(200+600)*0.512;
topoplot(squeeze(mean(mean2plot(:,start:endpoint),2)),chanlocs,'emarker2' ,{[12 48 49 19 32 56],'.','black'});
title(sprintf('Grand average of feedback trials (%s-%s ms)', num2str(300),num2str(600)))
caxis([-8,8])
colorbar()


