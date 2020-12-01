## Define functions for cohen's d
d <- function(t,n) {
  d <- t/(sqrt(n))
  return(d)}

se.d <-function(d,n){
  se<-sqrt((1/n)+((d^2)/(2*n)))
  names(se)<-"se for d"
  return(se)}

ci.d <- function(d, SE){
  ci.d <- c((d-1.96*SE),(d+1.96*SE))
  names(ci.d) <- c("95% CI lower", "95% CI upper")
  return(ci.d)
}


## define variables for info storage
lme_models <- list()
lme_null_models <- list()
BF_BIC <- list()
N.cases <- list()
N.controls <- list()
d_dx <- list()
ci_d <- list()
aic_main <- list(); 
bic_main <- list(); 
ll_main <- list()


aidata=aidata[aidata$DependentonPrimaryDrug=='1',]

covs = c("Age", "Sex", "Site","Niclast30days","DependentonPrimaryDrug")
sig_regions=c('AI_accumb')
## loop through sig regions
for (region in sig_regions)
{
  ## remove missings per data-subset (i.e. separately for each region)
  aidata.complete <- aidata[,c(covs,region)][which(complete.cases(aidata[,c(covs,region)])),]
  
  ## calculate number of cases and controls
  N.cases[[region]] <- as.numeric(table(aidata.complete$DependentonPrimaryDrug)["1"])
  N.controls[[region]] <- as.numeric(table(aidata.complete$DependentonPrimaryDrug)["0"])
  outcome <- aidata.complete[,region]
  lme_models[[region]] <- lme(outcome ~ Niclast30days + Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  lme_null_models[[region]] <- lme(outcome ~ Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  ## calculate cohen's d effect size
  d_dx[[region]] <- d(summary(lme_models[[region]])$tTable["Niclast30days","t-value"], N.cases[[region]])
  ## calcualte 95% CI around cohen's d
  ci_d[[region]] <- ci.d(d_dx[[region]], se.d(d_dx[[region]], N.cases[[region]]))
  
  ## extract model fit measures
  aic_main[[region]] <- AIC(lme_models[[region]])
  bic_main[[region]] <- BIC(lme_models[[region]])
  ll_main[[region]] <- lme_models[[region]]$logLik
  
  ## caculate BF via BIC approach
  BF_BIC[[region]] <- exp((BIC(lme_null_models[[region]]) - BIC(lme_models[[region]]))/2)
  
  ## clean
  rm(aidata.complete, outcome)
}
rm(region)

  #=====================================Extract model results==============================
data.models = lapply(lme_models, function(x){ summary(x)$tTable })
table <- data.frame(matrix(nrow=length(sig_regions[grep("AI", sig_regions)]), ncol=0))
table[,"region"] <- sub("AI_", "", sig_regions[grep("AI", sig_regions)])
table[,"Ntot"] <- (as.numeric(N.cases) + as.numeric(N.controls))
table[,"Ncases/controls"] <- paste(N.cases, "/", N.controls, sep="")
table[,"b_Niclast30days"] <- as.numeric(sapply(data.models, function(x){x[,"Value"]["Niclast30days"]}))
table[,"se_Niclast30days"] <- as.numeric(sapply(data.models, function(x){x[,"Std.Error"]["Niclast30days"]}))
table[,"t_Niclast30days"] <- as.numeric(sapply(data.models, function(x){x[,"t-value"]["Niclast30days"]}))
table[,"p_Niclast30days"] <- as.numeric(sapply(data.models, function(x){x[,"p-value"]["Niclast30days"]}))
table[,"cohensd.ci"] <- paste(round(as.numeric(d_dx),digits=3), "(", as.character(sapply(ci_d, function(p){ paste(round(p, digits=2), collapse=",", sep="")})), ")",sep="")
table[,"BF10"] <- as.numeric(BF_BIC)
savename <- 'Model_past30_nic.xls'
write.table(table, file=paste(working_dir,savename, sep=""), row.names = F, col.names=T, sep="\t")



