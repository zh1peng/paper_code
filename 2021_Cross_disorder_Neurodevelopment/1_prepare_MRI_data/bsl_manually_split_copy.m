addpath('/gpfs1/home/z/c/zcao4/matlab_tools/spm12')
subid=textread('bsl_log.txt','%s');
data_path='/gpfs2/scratch/acjulian/IMAGEN_files_still_need_BIDSifying/Task_BSL_Funcs/';
output_path='/gpfs2/scratch/zcao4/IMAGEN_MID_data/MID_BSL_EPI';
failed_subid={}
for indx=1:length(subid)
    subfolder=fullfile(data_path,subid{indx},'EPI_short_MID');
    output_subfolder=fullfile(output_path, subid{indx});
    if ~exist(output_subfolder,'dir')
        mkdir(output_subfolder)
    end
    try
        cd(output_subfolder)
        f2copy=spm_select('FPList', subfolder,'^wau.*\.nii.gz$');
        copyfile(f2copy, output_subfolder)
        f2unzip=spm_select('FPlist',output_subfolder,'^wau.*\.nii.gz$');
        gunzip(f2unzip)
        
        [filepath, fname, fexe]=fileparts(f2unzip);
        f2smooth=fname;
        P=spm_vol(f2smooth);
        Q=['s',f2smooth];
        s=[5,5,5];
        disp('smoothing......')
        spm_smooth(P,Q,s);
        
        f2split=Q;
        disp('splitting 4D ....')
        spm_file_split(f2split,output_subfolder);
        delete(f2split)
        delete(f2unzip)
        delete(f2smooth)
    catch
        failed_subid=[failed_subid, subid{indx}];
    end
end
