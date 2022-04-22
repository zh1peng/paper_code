
library(ape)
library(dplyr)
# data_path='F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data'
# lh_centroid_data=read.csv(sprintf('%s/lh_centroid.csv',data_path))
# centroid.l=lh_centroid_data %>%select(starts_with('centroid'))
# rownames(centroid.l)=lh_centroid_data$Row
# dist.mat=as.matrix(dist(centroid.l))
# w=1/dist.mat
# diag(w)= 0
# saveRDS(w, file = "F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data/Moran_w.rds")
w=readRDS('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data/Moran_w.rds')



### caculate moran's I for several brain maps
map.df=read.csv('F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data/brain_maps2test.csv')
map.df.l=map.df[1:34,-1]
map.df.moran=get_df_moran_idx(map.df.l)
write.csv(map.df.moran,sprintf('%s/Moran_I.csv',table_path))

#==========================================================================================
## create null brain new
## random brain but select I>0.01
## make sure it is a norm distribution around 0.03

set.seed(1990)
n_random_FPR=500000
random_brain_data=matrix(runif(34*n_random_FPR,min=-1,max=1),nrow=34)
random_brain_data=as.data.frame(random_brain_data)
colnames(random_brain_data)=paste0('null_brain_',c(1:n_random_FPR))

random_brain_moran=get_df_moran_idx(random_brain_data)
hist(brain_moran$moran_idx,breaks=100)

random_brain_moran.filtered=random_brain_moran[random_brain_moran$moran_idx>0,,drop=F]
random_brain_moran.filtered$prob=dnorm(random_brain_moran.filtered$moran_idx, mean = 0.03, sd = 0.01)


set.seed(1990)
sampled_moran=sample(random_brain_moran.filtered$moran_idx,
                     1000,
                     prob =random_brain_moran.filtered$prob)

# find their data

sampled_cols=rownames(random_brain_moran.filtered)[random_brain_moran.filtered$moran_idx %in% sampled_moran]
sampled_random_brain_data=random_brain_data[,sampled_cols]

# check moran idx
hist(get_df_moran_idx(sampled_random_brain_data)$moran_idx)

#rename colnames
colnames(sampled_random_brain_data)=paste0('null_',c(1:1000))
rownames(sampled_random_brain_data)=dk_labels.hem$label

# make some brain plots
p1= plot_brain(sampled_random_brain_data, 'null_1',show.legend = T)
p2= plot_brain(sampled_random_brain_data, 'null_2',show.legend = T)
p3= plot_brain(sampled_random_brain_data, 'null_79',show.legend = T)
p4= plot_brain(sampled_random_brain_data, 'null_107',show.legend = T)

grid.arrange(p1,p2,p3,p4,ncol=1)
saveRDS(sampled_random_brain_data, 'F:/Google Drive/post-doc/vitural_histology_revisit/clean_code/data/random_brain4FPR_uniform_subset1000.rds')




