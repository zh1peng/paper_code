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
null_type_level=c('random_gene',
                   'spin_brain')
null_type_label=c('Competitive null models',
                   'Self-contained null models')
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

#==1. analysis with simulation results
# Res_*csv file to analysis
analysis.path=file.path(result.path,sprintf('analysis_%s_%s_%s_%s_%s',atlas,rdonor,brain_type,gene_set_type,cor_type))
dir.create(analysis.path, showWarnings = F)


res.files=list(
  spin_brain=sprintf( '%s/Res_%s_%s_%s_%s_spin_brain_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type),
  random_gene=sprintf('%s/Res_%s_%s_%s_%s_random_gene_%s_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type,cor_type))
# read res.files
res.df.list=lapply(res.files, read.csv, stringsAsFactors = F)

#==1.1 Psig-G analysis
nest_by='geneSet'
pvals.nested=lapply(res.df.list, get_pvals_nested, nest_by=nest_by, heat_plot=F)
psig.list=lapply(pvals.nested, get_psig, if_fdr=F)
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
ggsave(file.path(analysis.path,'Psig_G_merged.png'),p12,width =15, height = 16 ,bg='white')

# co-expression analysis only makes sense when they are sorted by geneSet
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
p3=grid.arrange(grobs=coexp_res.plot.list$random_gene,
                ncol=2,
                top = textGrob("A. competitive null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))
p4=grid.arrange(grobs=coexp_res.plot.list$spin_brain,
                ncol=2,
                top = textGrob("B. Self-contained null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-G",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Co-expression",gp=gpar(fontsize=12,font=2)))
p34=grid.arrange(p3,p4,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))

ggsave(file.path(analysis.path,'Psig_G_coexp.png'),p34,width =10, height = 14 ,bg='white')


write.csv(coexp_res.report.list$spin_brain,file.path(analysis.path,'Psig_G_coexp_spin_brain.csv'),row.names = F)
write.csv(coexp_res.report.list$random_gene,file.path(analysis.path,'Psig_G_coexp_random_gene.csv'),row.names = F)


#==1.2 Psig-B analysis
nest_by='brain'
pvals.nested=lapply(res.df.list, get_pvals_nested, nest_by=nest_by, heat_plot=F)
psig.list=lapply(pvals.nested, get_psig, if_fdr=F)
p1=plot_violin_psig_list(psig.list = psig.list,
                         ylab2show = 'Psig-B',
                         title2show = 'A.',
                         title_adj = -0.086,
                         stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)
p2=plot_bar_psig_list(psig.list, 
                        ylab2show='Psig-B',
                        title2show = 'B.',
                        title_adj = -0.1,
                        stat_level = stat_level,
                         stat_label = stat_label,
                         null_type_level = null_type_level,
                         null_type_label = null_type_label)
p12=grid.arrange(p1,p2,ncol=1)
ggsave(file.path(analysis.path,'Psig_B_merged.png'),p12,width =15, height = 16 ,bg='white')
# get brain Moran's I and dip test stat of correlations with geneData
brain_info=get_brain_info(data_path='F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/data',
                          atlas=atlas,
                          rdonor=rdonor,
                          brain_type=brain_type,
                          method=cor_type)

brain_res.nested.list=lapply(psig.list, correlate_psig_with_info,info=brain_info,var2test='pos_neg_dist')
brain_res.report.list=lapply(brain_res.nested.list, report_res.nested)
brain_res.plot.list=lapply(brain_res.nested.list,plot_res.nested, 
                           annot_position=c(0.115,0.95),
                           xlim2show=c(-0.02,1),
                           ylim2show=c(-0.05,1))
brain_res.plot.list$random_gene=plot_res.nested(brain_res.nested.list$random_gene,
                                                annot_position=c(0.115,0.19),
                                                xlim2show=c(-0.02,1),
                                                ylim2show=c(-0.05,0.2))
p3=grid.arrange(grobs=brain_res.plot.list$random_gene,
                ncol=2,
                top = textGrob("A. competitive null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p4=grid.arrange(grobs=brain_res.plot.list$spin_brain,
                ncol=2,
                top = textGrob("B. Self-contained null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p34=grid.arrange(p3,p4,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))

ggsave(file.path(analysis.path,'Psig_B_BIdist.png'),p34,width =10, height = 14 ,bg='white')

write.csv(brain_res.report.list$spin_brain,file.path(analysis.path,'Psig_B_Bidist_spin_brain.csv'),row.names = F)
write.csv(brain_res.report.list$random_gene,file.path(analysis.path,'Psig_B_Bidist_random_gene.csv'),row.names = F)


# ============================modetest_stat================================
brain_res.nested.list=lapply(psig.list, correlate_psig_with_info,info=brain_info,var2test='modetest_stat')
brain_res.report.list=lapply(brain_res.nested.list, report_res.nested)
brain_res.plot.list=lapply(brain_res.nested.list,plot_res.nested,
                           annot_position=c(0,0.95),
                           xlim2show=c(-0.005,0.042),
                           ylim2show=c(-0.05,1))


p5=grid.arrange(grobs=brain_res.plot.list$random_gene,
                ncol=2,
                top = textGrob("A. competitive null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p6=grid.arrange(grobs=brain_res.plot.list$spin_brain,
                ncol=2,
                top = textGrob("B. Self-contained null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p56=grid.arrange(p5,p6,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))
ggsave(file.path(analysis.path,'Psig_B_BIdip.png'),p56,width =10, height = 14 ,bg='white')


write.csv(brain_res.report.list$spin_brain,file.path(analysis.path,'Psig_B_Bidip_spin_brain.csv'),row.names = F)
write.csv(brain_res.report.list$random_gene,file.path(analysis.path,'Psig_B_Bidip_random_gene.csv'),row.names = F)

# kval
brain_res.nested.list=lapply(psig.list, correlate_psig_with_info,info=brain_info,var2test='kval')
brain_res.report.list=lapply(brain_res.nested.list, report_res.nested)
brain_res.plot.list=lapply(brain_res.nested.list,plot_res.nested,
                           annot_position=c(0,1),
                           xlim2show=c(0,4),
                           ylim2show=c(-0.05,1))


p5=grid.arrange(grobs=brain_res.plot.list$random_gene,
                ncol=2,
                top = textGrob("A. competitive null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p6=grid.arrange(grobs=brain_res.plot.list$spin_brain,
                ncol=2,
                top = textGrob("B. Self-contained null models",gp=gpar(fontsize=16,font=1),x = -0.01, hjust = 0),
                left =textGrob("Psig-B",gp=gpar(fontsize=12,font=2),rot=90),
                bottom=textGrob("Bimodality",gp=gpar(fontsize=12,font=2)))
p56=grid.arrange(p5,p6,layout_matrix=rbind(matrix(rep(c(1,1,1),20),ncol=3),
                                           c(NA,NA,NA),
                                           matrix(rep(c(2,2,2),20),ncol=3)))
ggsave(file.path(analysis.path,'Psig_B_meanval.png'),p56,width =10, height = 14 ,bg='white')


write.csv(brain_res.report.list$spin_brain,file.path(analysis.path,'Psig_B_Bidip_spin_brain.csv'),row.names = F)
write.csv(brain_res.report.list$random_gene,file.path(analysis.path,'Psig_B_Bidip_random_gene.csv'),row.names = F)

