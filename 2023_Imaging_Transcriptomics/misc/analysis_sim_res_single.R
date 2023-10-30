project_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code'
source(sprintf('%s/functions/analysis_functions.R',project_path))
source(sprintf('%s/functions/data_functions.R',project_path))
source(sprintf('%s/functions/cor_functions.R',project_path))
sim_res_path=sprintf('%s/results',project_path)
result.path=sprintf('%s/reports',project_path)
dir.create(result.path, showWarnings = F)

atlas='desikan'
rdonor='r0.4'
brain_type='sim_spatial0.03'
gene_set_type='Sim'
cor_type='pearson'


null_type_level=c('random_gene_coexp',
                   'random_gene_subset',
                   'spin_random_mixed')
null_type_label=c('Null-coexp-gene',
                   'Null-brain-gene',
                   'Spin-random-mixed')
stat_level=c('mean',
            'meanabs',
            'meansqr',
            'maxmean',
            'median',
            'sig_n',
            'ks_orig',
            'ks_weighted')
stat_label=c('Mean',
            'Meanabs',
            'Meansqr',
            'Maxmean',
            'Median',
            'Sig Number',
            'KS',
            'Weighted KS')

# Res_*csv file to analysis
analysis.path=file.path(result.path,sprintf('analysis1_%s_%s_%s_%s_%s',atlas,rdonor,brain_type,gene_set_type,cor_type))
dir.create(analysis.path, showWarnings = F)


res.files=list(
  random_gene=sprintf( '%s/Res_%s_%s_%s_%s_random_gene_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type),
  random_gene_coexp=sprintf( '%s/Res_%s_%s_%s_%s_random_gene_coexp_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type),
  random_gene_subset=sprintf('%s/Res_%s_%s_%s_%s_random_gene_subset_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type),
  spin_random_mixed=sprintf('%s/Res_%s_%s_%s_%s_spin_random_mixed_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type))
# read res.files
res.df.list=lapply(res.files, read.csv, stringsAsFactors = F)

nest_by='geneSet'
pvals.nested=lapply(res.df.list, get_pvals_nested, nest_by=nest_by, heat_plot=F)
psig.list=lapply(pvals.nested, get_psig, if_fdr=F)

coexp_info=get_geneSetList_info(data_path=sprintf('%s/data',project_path),
                                 gs_type=gene_set_type,
                                 atlas=atlas,
                                 rdonor=rdonor)
coexp_res.nested.list=lapply(psig.list, correlate_psig_with_info,info=coexp_info,var2test='coexp_mean')
coexp_res.report.list=lapply(coexp_res.nested.list, report_res.nested)
coexp_res.plot.list=lapply(coexp_res.nested.list, 
                           plot_res.nested, 
                           xlim2show=c(-0.02,0.11),
                           annot_position=c(-0.01,0.5))
# random gene vs. random gene coexp matched
#===================================================================================================
null_type_level =c('random_gene','random_gene_coexp')
null_type_label =c('competitive Null models','Coexp-matched competitive Null models')
p1=plot_violin_psig_list(psig.list = psig.list,
                         ylab2show = 'Psig-G',
                         title2show = 'A.',
                         title_adj = -0.07,
                         stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)

p2=plot_bar_psig_list(psig.list, 
                        ylab2show='Psig-G',
                        title2show = 'B.',
                        title_adj = -0.1,
                        stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)

p12=grid.arrange(p1,p2,ncol=1)
ggsave(file.path(analysis.path,'Psig_G_Coexp-matched.png'),p12,width =15, height = 16 ,bg='white')


p3=grid.arrange(grobs=coexp_res.plot.list[[null_type_level[1]]],
                ncol=2,
                top = textGrob(sprintf("A. %s",null_type_label[1]),gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))
p4=grid.arrange(grobs=coexp_res.plot.list[[null_type_level[2]]],
                ncol=2,
                top = textGrob(sprintf("B. %s",null_type_label[2]),gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))


p34=grid.arrange(p3,p4,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))
ggsave(file.path(analysis.path,'Psig_G_coexp_coexp_matched.png'),p34,width =10, height = 14 ,bg='white')


# random gene vs. random gene subset
#===================================================================================================
null_type_level =c('random_gene','random_gene_subset')
null_type_label =c('competitive Null models','Brain competitive Null models')
p1=plot_violin_psig_list(psig.list = psig.list,
                         ylab2show = 'Psig-G',
                         title2show = 'A.',
                         title_adj = -0.07,
                         stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)

p2=plot_bar_psig_list(psig.list, 
                        ylab2show='Psig-G',
                        title2show = 'B.',
                        title_adj = -0.1,
                        stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)

p12=grid.arrange(p1,p2,ncol=1)
ggsave(file.path(analysis.path,'Psig_G_subset.png'),p12,width =15, height = 16 ,bg='white')




p3=grid.arrange(grobs=coexp_res.plot.list[[null_type_level[1]]],
                ncol=2,
                top = textGrob(sprintf("A. %s",null_type_label[1]),gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))
p4=grid.arrange(grobs=coexp_res.plot.list$random_gene_coexp,
                ncol=2,
                top = textGrob(sprintf("B. %s",null_type_label[2]),gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))


p34=grid.arrange(p3,p4,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))
ggsave(file.path(analysis.path,'Psig_G_coexp_subset.png'),p34,width =10, height = 14 ,bg='white')








