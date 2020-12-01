##============================== Dummy Code Substance =================================
library(plyr)
library(fastDummies)
aidata$Substance_name <- mapvalues(aidata$PrimaryDrug,from=c('0','1','2','3','4','5'), to=c('ctl','alc','nic','coc','met','can'))
aidata=dummy_cols(aidata, select_columns='Substance_name',remove_first_dummy = F)
LRdata$Substance_name <- mapvalues(LRdata$PrimaryDrug,from=c('0','1','2','3','4','5'), to=c('ctl','alc','nic','coc','met','can'))
LRdata=dummy_cols(LRdata, select_columns='Substance_name',remove_first_dummy = F)
covs=c("Age", "Sex","Substance_name", "Site",'Substance_name_ctl','Substance_name_alc','Substance_name_nic','Substance_name_coc','Substance_name_met','Substance_name_can')


#============================= Define functions for cohen's d================================
d <- function(df, t, n1, n2) {
  d <- t*(n1+n2)/(sqrt(n1*n2)*sqrt(df))
  return(d)}

se.d <-function(d,n1,n2){
  se<-sqrt((n1+n2)/(n1*n2)+(d^2)/(2*(n1+n2-2)))
  names(se)<-"se for d"
  return(se)}

ci.d <- function(d, SE){
  ci.d <- c((d-1.96*SE),(d+1.96*SE))
  names(ci.d) <- c("95% CI lower", "95% CI upper")
  return(ci.d)
}
##================================== run Model2 ========================================
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

## loop through all brain regions
for (region in colnames(aidata)[grep("AI", colnames(aidata))]){
  
  ## remove missings per data-subset (i.e. separately for each region)
  aidata.complete <- aidata[,c(covs,region)][which(complete.cases(aidata[,c(covs,region)])),]
  
  ## calculate number of cases and controls
  N.cases[[region]] <- as.numeric(table(aidata.complete$Substance_name)[substring(substance2test,nchar(substance2test)-2)])
  N.controls[[region]] <- as.numeric(table(aidata.complete$Substance_name)['ctl'])
  
  ## perform regression analysis on complete data (i.e. without missings)
  outcome <- aidata.complete[,region]
  lme_models[[region]] <- lme(outcome ~Substance_name_alc+Substance_name_nic+Substance_name_coc+Substance_name_met+Substance_name_can + Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  
  ## remove substance2test var to built the null model
  pred_var=c('Substance_name_alc','Substance_name_nic','Substance_name_coc','Substance_name_met','Substance_name_can',"Age", "Sex")
  remove_var=substance2test
  null_pred_var=setdiff(pred_var,remove_var)
  lme_null_models[[region]] <- update(lme_models[[region]],fixed.=as.formula(paste("outcome ~", paste(null_pred_var,collapse="+"))))
  
  
  ## calculate cohen's d effect size for diagnosis
  d_dx[[region]] <- d(summary(lme_models[[region]])$tTable[substance2test,"DF"], summary(lme_models[[region]])$tTable[substance2test,"t-value"], N.cases[[region]], N.controls[[region]])
  
  ## calcualte 95% CI around cohen's d
  ci_d[[region]] <- ci.d(d_dx[[region]], se.d(d_dx[[region]], N.cases[[region]], N.controls[[region]]))
  
  ## extract model fit measures
  aic_main[[region]] <- AIC(lme_models[[region]])
  bic_main[[region]] <- BIC(lme_models[[region]])
  ll_main[[region]] <- lme_models[[region]]$logLik
  BF_BIC[[region]] <- exp((BIC(lme_null_models[[region]]) - BIC(lme_models[[region]]))/2)
  
  ## clean
  rm(aidata.complete, outcome)
}
rm(region)




##========================================= Extract model results =================================================
data.models = lapply(lme_models, function(x){ summary(x)$tTable })
table <- data.frame(matrix(nrow=length(colnames(aidata)[grep("AI", colnames(aidata))]), ncol=0))

table[,"region"] <- sub("AI_", "", colnames(aidata)[grep("AI", colnames(aidata))])
table[,"Ntot"] <- (as.numeric(N.cases) + as.numeric(N.controls))
table[,"Ncases/controls"] <- paste(N.cases, "/", N.controls, sep="")
#table[,"b_int"] <- as.numeric(sapply(data.models, function(x){x[,"Value"]["(Intercept)"]}))
table[,"b_diag"] <- as.numeric(sapply(data.models, function(x){x[,"Value"][substance2test]}))
table[,"empty1"] <- NA
#table[,"se_int"] <- as.numeric(sapply(data.models, function(x){x[,"Std.Error"]["(Intercept)"]}))
table[,"se_diag"] <- as.numeric(sapply(data.models, function(x){x[,"Std.Error"][substance2test]}))
table[,"empty2"] <- NA
#table[,"t_int"] <- as.numeric(sapply(data.models, function(x){x[,"t-value"]["(Intercept)"]}))
table[,"t_diag"] <- as.numeric(sapply(data.models, function(x){x[,"t-value"][substance2test]}))

table[,"empty3"] <- NA
#table[,"p_int"] <- as.numeric(sapply(data.models, function(x){x[,"p-value"]["(Intercept)"]}))
table[,"p_diag"] <- as.numeric(sapply(data.models, function(x){x[,"p-value"][substance2test]}))

table[,"empty4"] <- NA
table[,"cohensd.ci"] <- paste(round(as.numeric(d_dx),digits=3), "(", as.character(sapply(ci_d, function(p){ paste(round(p, digits=2), collapse=",", sep="")})), ")",sep="")
table[,"BF10"] <- as.numeric(BF_BIC)
#	

## add column to indicate whether region survives mcc
mcc <- c(table[grep("thick|Thickness", table$region),][which(p.adjust(table[grep("thick|Thickness", table$region),]$p_diag, method="fdr") < 0.05),"region"],
         table[grep("surf|SurfArea", table$region),][which(p.adjust(table[grep("surf|SurfArea", table$region),]$p_diag, method="fdr") < 0.05),"region"],
         table[-grep("surf|thick|SurfArea|Thickness", table$region),][which(p.adjust(table[-grep("surf|thick|SurfArea|Thickness", table$region),]$p_diag, method="fdr") < 0.05),"region"])

cond <- paste("table$region=='",mcc, "'", collapse="|", sep="")
w <- eval(parse(text=paste("which(",cond,")", sep="")))
table$mcc <- "no"
table[w,"mcc"] <- "yes"


str2save=substring(substance2test,nchar(substance2test)-2)
## save table
write.table(table, file=paste('Model2_results_',str2save,'.xls',sep = ''), row.names = F, col.names=T, sep="\t")
rm(table)




