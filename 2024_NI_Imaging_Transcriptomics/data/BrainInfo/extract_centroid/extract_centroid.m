clear;clc
data_path='F:\Google Drive\post-doc\vitural_histology_revisit\revision_code\data\BrainInfo\extract_centroid'

atlas={'desikan','schaefer100','schaefer200','schaefer300'}
filename={'lh.aparc.annot',...
                'lh.Schaefer2018_100Parcels_7Networks_order.annot',...
                'lh.Schaefer2018_200Parcels_7Networks_order.annot',...
                'lh.Schaefer2018_300Parcels_7Networks_order.annot'}
 spere_file=fullfile(data_path, 'lh.sphere')
for atlas_i =1: length(atlas)
    atlas_label=atlas{atlas_i};
    annot_file=fullfile(data_path, filename{atlas_i});
    [label,centroid]=centroid_extraction_sphere(spere_file,annot_file);
    T = array2table(centroid, 'RowNames', label);
    
    % remove (first row) medial wall
    rowsToRemove = strcmp(T.Properties.RowNames, 'Background+FreeSurfer_Defined_Medial_Wall');
    T(rowsToRemove, :) = [];
    writetable(T, strcat(atlas{atlas_i},'_centroid.csv'),'WriteRowNames',1);
end
