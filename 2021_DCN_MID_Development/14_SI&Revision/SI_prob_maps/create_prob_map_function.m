function []=create_prob_map_function(data_path, output_path,subid, sub_cov,cov_names, con_num,sample_size,rep_time,emask)
%% addin emask
%% seperated by sex presuming sex is in 1st column of cov.
if ~exist(output_path)
    mkdir(output_path);
end
%msgbox('is the sex in the 1st columns in the cov?');
[male_cov,male_row]=select_group(sub_cov,2,'0',1); % column2 is the sex! here
male_id=subid(male_row);
[female_cov,female_row]=select_group(sub_cov,2,'1',1);
female_id=subid(female_row);

for n=1:rep_time
    %% 2nd level
    male_perm=randperm((length(male_id)));
    male_perm_indx=male_perm(1:sample_size);
    selected_male_id=male_id(male_perm_indx);
    selected_male_cov=male_cov(male_perm_indx,:);

    female_perm=randperm((length(female_id)));
    female_perm_indx=female_perm(1:sample_size);
    selected_female_id=female_id(female_perm_indx);
    selected_female_cov=female_cov(male_perm_indx,:);
    
    selected_id=[selected_male_id;selected_female_id];
    selected_cov=[selected_male_cov;selected_female_cov];
    output_1st=fullfile(output_path,num2str(n));
    level_2nd_emask(data_path, output_1st,selected_id, selected_cov,cov_names, con_num,emask)
%     level_2nd_emask(data_path, output_path,subid, sub_cov,cov_names, con_num,emask)


    clear selected_id selected_cov
    
    %% get corrected and uncorrected T value
    spmfile=fullfile(output_1st, 'SPM.mat')
    [un,fwe]=get_tvalues(spmfile); %1. make sure get_tvalues in the path 2. un==> uncorrected p=0.001||||||| fwe==> FWE corrected p=0.05;
    
    input=fullfile(output_1st, 'spmT_0001.nii');
    
    %%
    fwe_expression=sprintf('(i1>%0.2f) | (i1<-%0.2f)',fwe,fwe);
    fwe_output_name='binary_fwe';
    un_expression=sprintf('(i1>%0.2f) | (i1<-%0.2f)',un,un);
    un_output_name='binary_un';
    
    matlabbatch{1}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{1}.spm.util.imcalc.output = fwe_output_name;
    matlabbatch{1}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{1}.spm.util.imcalc.expression = fwe_expression;
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    matlabbatch{2}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{2}.spm.util.imcalc.output = un_output_name;
    matlabbatch{2}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{2}.spm.util.imcalc.expression = un_expression;
    matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
    matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
    %% fwe seperate
    fwe_pos_expression=sprintf('i1>%0.2f',fwe);
    fwe_pos_output_name='binary_fwe_pos';
    fwe_neg_expression=sprintf('i1<-%0.2f',fwe);
    fwe_neg_output_name='binary_fwe_neg';
    
    matlabbatch{3}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{3}.spm.util.imcalc.output = fwe_pos_output_name;
    matlabbatch{3}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{3}.spm.util.imcalc.expression = fwe_pos_expression;
    matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{3}.spm.util.imcalc.options.mask = 0;
    matlabbatch{3}.spm.util.imcalc.options.interp = 1;
    matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
    
    matlabbatch{4}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{4}.spm.util.imcalc.output = fwe_neg_output_name;
    matlabbatch{4}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{4}.spm.util.imcalc.expression = fwe_neg_expression;
    matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{4}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{4}.spm.util.imcalc.options.mask = 0;
    matlabbatch{4}.spm.util.imcalc.options.interp = 1;
    matlabbatch{4}.spm.util.imcalc.options.dtype = 4;
    
    %% uncorrected seperate
    un_pos_expression=sprintf('i1>%0.2f',un);
    un_pos_output_name='binary_un_pos';
    un_neg_expression=sprintf('i1<-%0.2f',un);
    un_neg_output_name='binary_un_neg';
    
    matlabbatch{5}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{5}.spm.util.imcalc.output = un_pos_output_name;
    matlabbatch{5}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{5}.spm.util.imcalc.expression = un_pos_expression;
    matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{5}.spm.util.imcalc.options.mask = 0;
    matlabbatch{5}.spm.util.imcalc.options.interp = 1;
    matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
    
    matlabbatch{6}.spm.util.imcalc.input = cellstr(input);
    matlabbatch{6}.spm.util.imcalc.output = un_neg_output_name;
    matlabbatch{6}.spm.util.imcalc.outdir = cellstr(output_1st);
    matlabbatch{6}.spm.util.imcalc.expression = un_neg_expression;
    matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{6}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{6}.spm.util.imcalc.options.mask = 0;
    matlabbatch{6}.spm.util.imcalc.options.interp = 1;
    matlabbatch{6}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch);
    clear matlabbatch
    
    
    
    
    
    all_fwe_binary(n,1)=cellstr(fullfile(output_1st,[fwe_output_name,'.nii']));
    all_fwe_pos_binary(n,1)=cellstr(fullfile(output_1st,[fwe_pos_output_name,'.nii']));
    all_fwe_neg_binary(n,1)=cellstr(fullfile(output_1st,[fwe_neg_output_name,'.nii']));
    all_un_binary(n,1)=cellstr(fullfile(output_1st,[un_output_name,'.nii']));
    all_un_pos_binary(n,1)=cellstr(fullfile(output_1st,[un_pos_output_name,'.nii']));
    all_un_neg_binary(n,1)=cellstr(fullfile(output_1st,[un_neg_output_name,'.nii']));
    
end

%% final results
matlabbatch{1}.spm.util.imcalc.input = all_fwe_binary;
matlabbatch{1}.spm.util.imcalc.output = 'all_fwe_binary';
matlabbatch{1}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.util.imcalc.input = all_fwe_pos_binary;
matlabbatch{2}.spm.util.imcalc.output = 'all_fwe_pos_binary';
matlabbatch{2}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{2}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
matlabbatch{3}.spm.util.imcalc.input = all_fwe_neg_binary;
matlabbatch{3}.spm.util.imcalc.output = 'all_fwe_neg_binary';
matlabbatch{3}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{3}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 4;
matlabbatch{4}.spm.util.imcalc.input = all_un_binary;
matlabbatch{4}.spm.util.imcalc.output = 'all_un_binary';
matlabbatch{4}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{4}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{4}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{4}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{4}.spm.util.imcalc.options.mask = 0;
matlabbatch{4}.spm.util.imcalc.options.interp = 1;
matlabbatch{4}.spm.util.imcalc.options.dtype = 4;
matlabbatch{5}.spm.util.imcalc.input =  all_un_pos_binary;
matlabbatch{5}.spm.util.imcalc.output = 'all_un_pos_binary';
matlabbatch{5}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{5}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{5}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{5}.spm.util.imcalc.options.mask = 0;
matlabbatch{5}.spm.util.imcalc.options.interp = 1;
matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
matlabbatch{6}.spm.util.imcalc.input = all_un_neg_binary;
matlabbatch{6}.spm.util.imcalc.output = 'all_un_neg_binary';
matlabbatch{6}.spm.util.imcalc.outdir = cellstr(output_path);
matlabbatch{6}.spm.util.imcalc.expression = 'mean(X)';
matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{6}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{6}.spm.util.imcalc.options.mask = 0;
matlabbatch{6}.spm.util.imcalc.options.interp = 1;
matlabbatch{6}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch
end



   function [u,uc]=get_tvalues(spm_dir)
%% calculate t value for a given SPM.mat file 
% u is for uncorrected 0.001
% u1 is for fwe corrected 0.05
load(spm_dir)

df   = [1 SPM.xX.erdf];
STAT = 'T';
R    = SPM.xVol.R;
S    = SPM.xVol.S;
n=1;
[u] = spm_u(0.001,df,STAT);
[uc] = spm_uc(0.05,df,STAT,R,n,S);

end 






function []=level_2nd_emask(data_path, output_path,subid, sub_cov,cov_names, con_num,emask)
%Updated 20/06/2016 emask added in
spm_jobman('initcfg')
spm('Defaults','fMRI')
spm_get_defaults('cmdline',true)

if ~exist(output_path)
    mkdir(output_path)
end
if con_num>9
   sub_cons=fullfile(data_path,subid,['con_00',num2str(con_num),'.nii']);
else
    sub_cons=fullfile(data_path,subid,['con_000',num2str(con_num),'.nii']);
end


matlabbatch=[];
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(output_path);
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans=sub_cons;
    for k=1:length(cov_names)
        matlabbatch{1}.spm.stats.factorial_design.cov(k).cname = cov_names{k};
        matlabbatch{1}.spm.stats.factorial_design.cov(k).c =sub_cov(:,k);
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCC = 1;
    end
matlabbatch{1}.spm.stats.factorial_design.masking.im=1;
if emask~=0
matlabbatch{1}.spm.stats.factorial_design.masking.em =cellstr(emask);
end
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'pos';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
spm_jobman('run',matlabbatch);
matlabbatch=[];
end

