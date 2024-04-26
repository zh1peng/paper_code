# report odd simulated brain maps

for (stat2test in c('meanabs','meansqr','maxmean','mean','median','ks_orig','ks_weighted')){
odd.brain=psig.list$spin_brain %>% filter(.data[[stat2test]]>0.9)
#print(odd.brain$brain_id)
print(paste(stat2test, length(odd.brain$brain_id),sep=': '))
}



# find null_brain showing all significant correlation with random genes
 odd.brain=psig.list$spin_brain %>% filter(meanabs>.9)
 print(odd.brain$brain_id)
 odd_idx=as.numeric(stringr::str_remove(odd.brain$brain_id,'brain_'))
 odd.brain=load_BrainDat(atlas=atlas,type=brain_type,col_idx = 21)
 
 gene_data=load_GeneExp(atlas=atlas,rdonor = rdonor)
 cor2plot=corr_brain_gene(gene_data,odd.brain,method='pls1')
 cor2plot2=corr_brain_gene(gene_data,odd.brain,method='pearson')
 cor(cor2plot,cor2plot2)
 
 hist(cor2plot)

kurtval=moments::kurtosis(cor2plot)

perm.id=load_permid()
odd.brain.null=generate_null_brain_data(odd.brain,perm.id)
cor2plot.null=corr_brain_gene(gene_data,odd.brain.null[,c(1:200)],method='pls1w')

moments::kurtosis(cor2plot.null)

hist(cor2plot.null)



# Generate dataset1 (Normal distribution)
dataset1 <- rnorm(10000, mean = 0, sd = 1)
# Generate dataset2 (Bimodal distribution)
dataset2 <- c(rnorm(5000, mean = -2, sd = 1), rnorm(5000, mean = 2, sd = 1))
# Calculate skewness
library(moments)
skewness_dataset1 <- skewness(dataset1)
skewness_dataset2 <- skewness(dataset2)
# Calculate kurtosis
kurtosis_dataset1 <- kurtosis(dataset1)
kurtosis_dataset2 <- kurtosis(dataset2)
# Print the results
print(paste("Dataset1 - Skewness:", skewness_dataset1, "Kurtosis:", kurtosis_dataset1))
print(paste("Dataset2 - Skewness:", skewness_dataset2, "Kurtosis:", kurtosis_dataset2))



set.seed(101)
X <- replicate(2, rnorm(100))
y <- 0.6*X[,1] + 0.7*X[,2] + rnorm(100)
X <- apply(X, 2, scale)
y <- scale(y)


# NIPALS (PLS1)
u <- crossprod(X, y)
u <- u/drop(sqrt(crossprod(u)))         # X weights
t  <- X%*%u
p <- crossprod(X, t)/drop(crossprod(t)) # X loadings

coef <- u*drop(coef(lm(y~0+t))) # coefficient


coef(lm(y~0+X)) 
cor(cor(X,y), coef(plsr.model))
cor(u, coef(plsr.model))

# using pls package
library(pls)
plsr.model=plsr(y ~ X, ncomp=1)
coef(plsr.model)
loading.weights(plsr.model)
# colMeans(df.gi)
# df.li=get_column_localmoran(odd.brain.df)
# mean(apply(df.li, 2, sd))
# 
# #res.df.list$random_brain %>% filter(brain_id ==174) %>% select(ends_with('meanabs'))
# res.df.list$random_gene %>% filter(null_brain ==174) %>% select(ends_with('meanabs'))
# res.df$random_brain %>% filter(brain_id ==174) %>% select(ends_with('meanabs'))
# 
# gene_data=load_dk_gene_data(type = 'abagen_r04')
# brain_data=get_brain_data(type='random_spatial003',col_idx =970)
# brain_data=get_brain_data(type='random_spatial003',col_idx =1)
# brain_data=get_brain_data(type='real_brain',col_idx =5)
# 
# cor2plot=corr_brain_gene(gene_data,brain_data)
# hist(cor2plot)
# 
# A=multimode::modetest(cor2plot,method='HH')
# multimode::modetest(cor2plot,method='SI')
# multimode::modetest(cor2plot,method='ACR')
# diptest::dip(cor2plot)
# get_pos_neg_dist(cor2plot)
# 
# oddbrain=get_brain_data(type='random_spatial003',col_idx = odd.brain$null_brain)
# pc1=get_brain_data(type='real_brain',col_idx = 7)
# cor(pc1, oddbrain)