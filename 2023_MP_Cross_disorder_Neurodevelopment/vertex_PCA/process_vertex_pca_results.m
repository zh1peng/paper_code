function []=process_vertex_pca_results(mat_dir, mat_name)
mat2load=fullfile(mat_dir,[mat_name,'.mat'])
load(mat2load)
mkdir(fullfile(mat_dir,mat_name))
cd (fullfile(mat_dir,mat_name))

addpath('F:\Google Drive\post-doc\p_factor\IMAGEN_FU2_vertex')
addpath(genpath('F:\Google Drive\post-doc\p_factor\analysis_codes\surfstat'))
s = SurfStatReadSurf( { 'F:\Google Drive\post-doc\p_factor\IMAGEN_FU2_vertex\fsaverage\lh.pial',...
                                    'F:\Google Drive\post-doc\p_factor\IMAGEN_FU2_vertex\fsaverage\rh.pial'} );
SurfStatView( pc1_data', s, 'PC1' );
colormap(flipud(autumn))
saveas(gcf,'PC1.png')
SurfStatView( pc2_data', s, 'PC2' );
colormap jet
SurfStatColLim( [-0.5, 0.5] )
saveas(gcf,'PC2.png')

SurfStatView(pc3_data', s, 'PC3' );
colormap jet
SurfStatColLim( [-0.5, 0.5] )
saveas(gcf,'PC3.png')

% average loading within atlas
pc1_T=my_surf2parcel(pc1_data);
writetable(pc1_T,'pc1_vertex.csv','WriteRowNames',1)
pc2_T=my_surf2parcel(pc2_data);
writetable(pc2_T,'pc2_vertex.csv','WriteRowNames',1)
pc3_T=my_surf2parcel(pc3_data);
writetable(pc3_T,'pc3_vertex.csv','WriteRowNames',1)


% plot VAF

pca_variance=pca_sdev.^2
pca_vaf=pca_variance/sum(pca_variance)*100;
figure
plot(pca_vaf(1:10))
xlabel('PCs')
ylabel('VAF (%)')
saveas(gcf,'PCA_VAF.png')
end


