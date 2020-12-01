lme_models <- list()
resid_data <- list()
LRdata.complete = LRdata[complete.cases(LRdata[,grep("Age|Sex|Site|AI_", colnames(LRdata))]),]
for (region in colnames(LRdata.complete)[grep("L_|R_", colnames(LRdata.complete))] ){
  outcome <- LRdata.complete[,region]
  lme_models[[region]] <- lme(outcome ~ Sex + Age, random=~1|Site, data = LRdata.complete, method="ML", na.action="na.fail")
  resid_data[[region]] <- resid(lme_models[[region]], type='normalized')
  rm(outcome)
}
rm(region)
resid_LRdata=cbind(LRdata.complete[,c("Age", "Sex", "Site","Dependentanydrug","PrimaryDrug")],data.frame(resid_data))



