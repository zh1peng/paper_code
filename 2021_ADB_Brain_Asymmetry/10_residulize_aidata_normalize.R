lme_models <- list()
resid_data <- list()
aidata.complete = aidata[complete.cases(aidata[,grep("Age|Sex|Site|AI_", colnames(aidata))]),]
for (region in colnames(aidata.complete)[grep("AI", colnames(aidata.complete))] ){
  outcome <- aidata.complete[,region]
  lme_models[[region]] <- lme(outcome ~ Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  resid_data[[region]] <- resid(lme_models[[region]], type="Normalize")
  rm(outcome)
}
rm(region)
resid_aidata=cbind(aidata.complete[,c("Age", "Sex", "Site","Dependentanydrug","PrimaryDrug")],data.frame(resid_data))



