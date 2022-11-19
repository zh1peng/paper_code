results_path='F:\Google Drive\post-doc\p_factor\IMAGEN_BSL_vertex\vertex_pca_results'
addpath(results_path)
for n=[0,5,10,15, 20,25]
mat_name=sprintf('pca_ct_fwhm%d_resid_lm',n);
process_vertex_pca_results(results_path, mat_name)
end