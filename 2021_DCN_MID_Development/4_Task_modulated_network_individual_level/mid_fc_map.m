function con_PPI_final=mid_fc_map(SPM, mnis,name2save,sub_tag)
% just work for MID task
cd(SPM.swd)
my_con_nr=2;
my_radius=5;
my_name='ROIs';
% Use original mask order and sort later
for mni_i=1:size(mnis,1)
 disp(['Sub:', num2str(sub_tag), '~~~~ROI: ', num2str(mni_i)]);
my_xyz=mnis(mni_i,:);
[~,xY]=my_spm_regions(SPM, my_xyz, my_name, my_con_nr, my_radius);
PPI_no=simple_PPI(SPM, xY,[1,1,1],'no_win');
PPI_big=simple_PPI(SPM, xY,[3,1,1],'big_win');
% PPI_P_no(:,mni_i)=PPI_no.P;
% PPI_Y_no(:,mni_i)=PPI_no.Y;
% PPI_ppi_no(:,mni_i)=PPI_no.ppi;
% PPI_P_big(:,mni_i)=PPI_big.P;
% PPI_Y_big(:,mni_i)=PPI_big.Y;
% PPI_ppi_big(:,mni_i)=PPI_big.ppi;
PPI_matrix(:,:,mni_i)=[PPI_no.P,PPI_big.P,PPI_no.Y,PPI_no.ppi,PPI_big.ppi,ones(size(PPI_no.Y))]; % PPI_big.Y==PPI.no.Y
Y_ts(:,mni_i)=PPI_no.Y;
PPI_no=[];
PPI_big=[];
xY=[];
end

for roii=1:size(mnis,1)
    for roij=1:size(mnis,1)
        if roij==roii
                beta_PPI(:,roii,roij) = NaN(6,1);
                con_PPI(roii,roij) = NaN;
        else
                [b,~] = regress(zscore(Y_ts(:,roij)),squeeze(PPI_matrix(:,:,roii)));
                beta_PPI(:,roii,roij) = b;
                con_PPI(roii,roij) = [0 0 0 -1 1 0]*b;
        end
    end
end

% A=[1 2 3;4 5 6;7 8 9]
% rot90(fliplr(A),1)
con_PPI_flip=rot90(fliplr(con_PPI),1);
con_PPI_final=bsxfun(@plus, con_PPI_flip, con_PPI)./2;
% con_PPI_final(isnan(con_PPI_final))=1;
str2save=fullfile(SPM.swd,name2save);
save(str2save, 'con_PPI_final');
other_mat=fullfile(SPM.swd,'PPI_fc_all_info.mat');
save(other_mat,'con_PPI','beta_PPI','Y_ts','PPI_matrix')
end