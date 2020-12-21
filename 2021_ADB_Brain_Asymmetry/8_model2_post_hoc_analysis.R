library(nlme)
library(emmeans)
covs=c("Age", "Sex","Substance_name", "Site")
for (region in c('AI_accumb')){
  
  ## remove missings per data-subset (i.e. separately for each region)
  aidata.complete <- aidata[,c(covs,region)][which(complete.cases(aidata[,c(covs,region)])),]
  outcome <- aidata.complete[,region]
  lme_model <- lme(outcome ~Substance_name + Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  tmp=emmeans(lme_model, list(pairwise ~ Substance_name), adjust = "fdr")
  result_df=data.frame(tmp$`pairwise differences of Substance_name`)
  df2save=result_df %>% select(contrast,t.ratio,p.value)
  df2save$t.ratio=round(df2save$t.ratio,3)
  df2save$p.value=round(df2save$p.value,3)
write.csv(df2save,paste('post_hoc_',region,'.csv',sep=''))}



