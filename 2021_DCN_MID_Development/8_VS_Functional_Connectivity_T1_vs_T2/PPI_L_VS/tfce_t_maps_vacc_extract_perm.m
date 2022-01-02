cd('/gpfs1/home/z/c/zcao4/IMAGEN_MID_analysis/PPI_LME_new_subid_unbalanced/PPI_L_VS/perm_time')
addpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12')
nii_info=dir('int_t_vals*.nii')
sample_info=spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/sample_nii.nii')
sample_info.descrip='sample_spm_file'
mask_img=spm_read_vols(spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/matched_GM_mask.nii'));
pos_max=[]
neg_max=[]
for i=1:length(nii_info)
    disp(num2str(i))
    nii2read=nii_info(i).name;
    [~,nii_name,~]=fileparts(nii2read);
  
    vol_info=spm_vol(nii2read);
  img=spm_read_vols(vol_info);
  tmp_img=img.*mask_img;
%   pos
  tmp_img(tmp_img<0)=0;
 pos_tfced = matlab_tfce_transform(tmp_img,2,0.5,18,0.1);
pos_max=[pos_max,max(pos_tfced,[],'all')];
% sample_info.fname=strcat(nii_name,'_tfced_pos.nii');
 %spm_write_vol(sample_info,pos_tfced);
 
%   neg 
  tmp_img=img.*mask_img;
  tmp_img(tmp_img>0)=0;
  neg_tfced = matlab_tfce_transform(abs(tmp_img),2,0.5,18,0.1);
%  sample_info.fname=strcat(nii_name,'_tfced_neg.nii');
  %spm_write_vol(sample_info,neg_tfced);
neg_max=[neg_max, max(neg_tfced,[],'all')];
end
pos_max=sort(pos_max)
neg_max=sort(neg_max)
save('int_perm_result.mat','pos_max','neg_max')



clear 
nii_info=dir('t1_vals*.nii')
sample_info=spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/sample_nii.nii')
sample_info.descrip='sample_spm_file'
mask_img=spm_read_vols(spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/matched_GM_mask.nii'));
pos_max=[]
neg_max=[]
for i=1:length(nii_info)
    disp(num2str(i))
    nii2read=nii_info(i).name;
    [~,nii_name,~]=fileparts(nii2read);
  
    vol_info=spm_vol(nii2read);
  img=spm_read_vols(vol_info);
  tmp_img=img.*mask_img;
%   pos
  tmp_img(tmp_img<0)=0;
 pos_tfced = matlab_tfce_transform(tmp_img,2,0.5,18,0.1);
pos_max=[pos_max,max(pos_tfced,[],'all')];
% sample_info.fname=strcat(nii_name,'_tfced_pos.nii');
 %spm_write_vol(sample_info,pos_tfced);
 
%   neg 
  tmp_img=img.*mask_img;
  tmp_img(tmp_img>0)=0;
  neg_tfced = matlab_tfce_transform(abs(tmp_img),2,0.5,18,0.1);
%  sample_info.fname=strcat(nii_name,'_tfced_neg.nii');
  %spm_write_vol(sample_info,neg_tfced);
neg_max=[neg_max, max(neg_tfced,[],'all')];
end
pos_max=sort(pos_max)
neg_max=sort(neg_max)
save('t1_perm_result.mat','pos_max','neg_max')




clear pos_max neg_max

nii_info=dir('t2_vals*.nii')
sample_info=spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/sample_nii.nii')
sample_info.descrip='sample_spm_file'
mask_img=spm_read_vols(spm_vol('/gpfs1/home/z/c/zcao4/matlab_tools/spm12/matched_GM_mask.nii'));
pos_max=[]
neg_max=[]
for i=1:length(nii_info)
    disp(num2str(i))
    nii2read=nii_info(i).name;
    [~,nii_name,~]=fileparts(nii2read);
  
    vol_info=spm_vol(nii2read);
  img=spm_read_vols(vol_info);
  tmp_img=img.*mask_img;
%   pos
  tmp_img(tmp_img<0)=0;
 pos_tfced = matlab_tfce_transform(tmp_img,2,0.5,18,0.1);
pos_max=[pos_max,max(pos_tfced,[],'all')];
% sample_info.fname=strcat(nii_name,'_tfced_pos.nii');
 %spm_write_vol(sample_info,pos_tfced);
 
%   neg 
  tmp_img=img.*mask_img;
  tmp_img(tmp_img>0)=0;
  neg_tfced = matlab_tfce_transform(abs(tmp_img),2,0.5,18,0.1);
%  sample_info.fname=strcat(nii_name,'_tfced_neg.nii');
  %spm_write_vol(sample_info,neg_tfced);
neg_max=[neg_max, max(neg_tfced,[],'all')];
end
pos_max=sort(pos_max)
neg_max=sort(neg_max)
save('t2_perm_result.mat','pos_max','neg_max')








