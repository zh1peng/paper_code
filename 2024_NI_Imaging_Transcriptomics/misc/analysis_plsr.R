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

analysis.path=file.path(result.path,sprintf('pslr_%s_%s_%s_%s_%s',atlas,rdonor,brain_type,gene_set_type,cor_type))
dir.create(analysis.path, showWarnings = F)


res.files=list(
  spin_brain=sprintf( '%s/Res_%s_%s_%s_%s_spin_brain_pearson_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type),
  random_gene=sprintf('%s/Res_%s_%s_%s_%s_random_gene_pearson_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type),
  plsr_spin_brain=sprintf('%s/Res_%s_%s_%s_%s_spin_brain_PLSR_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type),
  plsr_random_gene=sprintf('%s/Res_%s_%s_%s_%s_random_gene_PLSR_sim1000.csv',sim_res_path,atlas,rdonor,brain_type,gene_set_type))
# read res.files
res.df.list=lapply(res.files, read.csv, stringsAsFactors = F)

nest_by='geneSet'
pvals.nested=lapply(res.df.list, get_pvals_nested, nest_by=nest_by, heat_plot=F)
psig.list=lapply(pvals.nested, get_psig, if_fdr=F)

psig.list$random_gene=merge(psig.list$random_gene, psig.list$plsr_random_gene, by=c('geneSet'))
psig.list$spin_brain=merge(psig.list$spin_brain, psig.list$plsr_spin_brain, by=c('geneSet'))
psig.list$plsr_random_gene=NULL
psig.list$plsr_spin_brain=NULL

null_type_level=c('random_gene',
                   'spin_brain')
null_type_label=c('competitive Null Models',
                   'Self-Contained Null Models')
stat_level=c('mean',
            'meanabs',
            'meansqr',
            'maxmean',
            'median',
            'sig_n',
            'ks_orig',
            'ks_weighted',
            'PLSR')
stat_label=c('Mean',
            'Meanabs',
            'Meansqr',
            'Maxmean',
            'Median',
            'Sig Number',
            'KS',
            'Weighted KS',
            'PLSR')


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



nest_by='brain'
pvals.nested=lapply(res.df.list, get_pvals_nested, nest_by=nest_by, heat_plot=F)
psig.list=lapply(pvals.nested, get_psig, if_fdr=F)
psig.list$random_gene=merge(psig.list$random_gene, psig.list$plsr_random_gene, by=c('brain_id'))
psig.list$spin_brain=merge(psig.list$spin_brain, psig.list$plsr_spin_brain, by=c('brain_id'))
psig.list$plsr_random_gene=NULL
psig.list$plsr_spin_brain=NULL
p1=plot_violin_psig_list(psig.list = psig.list,
                         ylab2show = 'Psig-B',
                         title2show = 'A.',
                         title_adj = -0.07,
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
