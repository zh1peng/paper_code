setwd('/gpfs1/home/z/c/zcao4/PCA_project/IMAGEN_FU2_vertex')
library('freesurferformats')
library('nlme')
library('dplyr')
library('R.matlab')
library('fastDummies')

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Need to specify fwhm (0, 5, 25) and if_resid(none/lme/lm)", call. = FALSE)
}

fwhm = args[1]
if_resid = args[2]

print('=============reading vertex data===========')
lh_ct = read.fs.mgh(sprintf('lh.all.thickness.fwhm%s.mgh', fwhm))
vx_n = dim(lh_ct)[1]
lh_ct_data = t(drop(lh_ct))
rh_ct = read.fs.mgh(sprintf('rh.all.thickness.fwhm%s.mgh', fwhm))
rh_ct_data = t(drop(rh_ct))

lh_ct_mask = read.fs.label('lh.cortex.label')
rh_ct_mask = read.fs.label('rh.cortex.label')
ct_data = cbind(lh_ct_data[, lh_ct_mask], rh_ct_data[, rh_ct_mask])
rm(lh_ct, lh_ct_data, rh_ct, rh_ct_data)


##===========================Set numeric cols==============================================
print('=============reading covs===========')
df = read.csv('cov_info.csv')
#set data format
table(apply(df, 2, class)) ## all columns are of class character

df2 <- df
for (col in colnames(df2)[grep("Age|ICV|mod_PDS", colnames(df))]) {
  df2[, col] <- as.character(df[, col])
  df2[, col] <- as.numeric(df2[, col])
}
## checked, values don't change after setting to numeric. FALSE means no changes
TRUE %in% apply(df == df2, 2, function(x) {
  FALSE %in% x
})
df <- df2
rm(col, df2)

df$ICV = df$ICV / 1000000 #rescale ICV

df$Sex = as.factor(as.character(df$Sex))
df$Sex = relevel(df$Sex, '0')

df$Site = as.factor(as.character(df$Site))
fac_num = df %>%
  group_by(Site) %>%
  dplyr::summarise(n_rows = length(Site)) %>%
  arrange(n_rows) %>%
  mutate_if(is.factor, as.character) # as.charater(factor) will give only numbers!!! as.character(fac_num[dim(fac_num)[1],1])
# set the largest group as reference
df$Site = relevel(df$Site, as.character(fac_num[dim(fac_num)[1], 1]))

if (if_resid == 'lme') {
  print('=============residulizing with LME===========')
  resid_data = list()
  for (vx_i in c(1:dim(ct_data)[2])) {
    region = sprintf('vx%s', vx_i)
    df$outcome <- ct_data[, vx_i]
    lme_model <- lme(
      outcome ~ Sex + Age + ICV,
      random =  ~ 1 | Site,
      data = df,
      control = lmeControl(opt = "optim")
    )
    resid_data[[region]] <- resid(lme_model, type = 'response')
  }
  ct_data.resid = as.matrix(data.frame(resid_data))
  print('==============runing PCA with LME resid===========')
  pca = prcomp(ct_data.resid, center = T, scale. = T)
} else if (if_resid == 'lm') {
  print('=============residulizing with LM===========')
  X = df %>% select(Age, Sex, ICV, Site)
  X = dummy_cols(X,
                select_columns = 'Site',
                remove_first_dummy = TRUE,
                remove_selected_columns = TRUE
  )
  X = data.matrix(X)
  y = ct_data
  A <- solve(t(X) %*% X)
  beta <- A %*% t(X) %*% y
  ct_data.resid <- y - X %*% beta
  print('==============runing PCA with LM resid===========')
  pca = prcomp(ct_data.resid, center = T, scale. = T)
} else if (if_resid == 'none') {
  print('==============runing PCA without resid===========')
  pca = prcomp(ct_data, center = T, scale. = T)
}
#

pca$rotation = pca$rotation %*% diag(as.numeric(diag(pca$rotation) > 0) *
                                       2 - 1)

print('============saving mat files=============')


lr_data_idx = c(lh_ct_mask, rh_ct_mask + vx_n)
#PC1
pc1.loadings = pca$rotation[, 1] * pca$sdev[1]
pc1.data = numeric(vx_n * 2)
for (idx in c(1:length(lr_data_idx))) {
  pc1.data[lr_data_idx[idx]] <- pc1.loadings[idx]
}



# PC2
pc2.loadings = pca$rotation[, 2] * pca$sdev[2]
pc2.data = numeric(vx_n * 2)
for (idx in c(1:length(lr_data_idx))) {
  pc2.data[lr_data_idx[idx]] <- pc2.loadings[idx]
}


# PC3
pc3.loadings = pca$rotation[, 3] * pca$sdev[3]
pc3.data = numeric(vx_n * 2)
for (idx in c(1:length(lr_data_idx))) {
  pc3.data[lr_data_idx[idx]] <- pc3.loadings[idx]
}

pca_sdev = pca$sdev

writeMat(
  sprintf("pca_ct_fwhm%s_resid_%s.mat", fwhm, if_resid),
  pc1_data = pc1.data,
  pc2_data = pc2.data,
  pc3_data = pc3.data,
  pca_sdev = pca_sdev
)
