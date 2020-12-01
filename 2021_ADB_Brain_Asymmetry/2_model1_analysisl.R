##=============================
# Model1: dependent vs. non depedent

library(nlme)
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


#================================================LME model============================================
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
cov = c("Age", "Sex", "Dependentanydrug", "Site")

## loop through all brain regions
for (region in colnames(aidata)[grep("AI", colnames(aidata))] )
{
  ## remove missings per data-subset (i.e. separately for each region)
  aidata.complete <- aidata[,c(cov,region)][which(complete.cases(aidata[,c(cov,region)])),]
  
  ## calculate number of cases and controls
  N.cases[[region]] <- as.numeric(table(aidata.complete$Dependentanydrug)["1"])
  N.controls[[region]] <- as.numeric(table(aidata.complete$Dependentanydrug)["0"])
  
  ## perform regression analysis on complete data (i.e. without missings)
  outcome <- aidata.complete[,region]
  lme_models[[region]] <- lme(outcome ~ Dependentanydrug + Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  lme_null_models[[region]] <- lme(outcome ~ Sex + Age, random=~1|Site, data = aidata.complete, method="ML", na.action="na.fail")
  
  ## calculate cohen's d effect size for diagnosis
  d_dx[[region]] <- d(summary(lme_models[[region]])$tTable["Dependentanydrug1","DF"], 
                      summary(lme_models[[region]])$tTable["Dependentanydrug1","t-value"], 
                      N.cases[[region]], 
                      N.controls[[region]])
  
  ## calcualte 95% CI around cohen's d
  ci_d[[region]] <- ci.d(d_dx[[region]], se.d(d_dx[[region]], N.cases[[region]], N.controls[[region]]))
  
  ## extract model fit measures
  aic_main[[region]] <- AIC(lme_models[[region]])
  bic_main[[region]] <- BIC(lme_models[[region]])

  ## caculate BF via BIC approach
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
table[,"b_diag"] <- as.numeric(sapply(data.models, function(x){x[,"Value"]["Dependentanydrug1"]}))
table[,"empty1"] <- NA
table[,"se_diag"] <- as.numeric(sapply(data.models, function(x){x[,"Std.Error"]["Dependentanydrug1"]}))
table[,"empty2"] <- NA
table[,"t_diag"] <- as.numeric(sapply(data.models, function(x){x[,"t-value"]["Dependentanydrug1"]}))
table[,"empty3"] <- NA
table[,"p_diag"] <- as.numeric(sapply(data.models, function(x){x[,"p-value"]["Dependentanydrug1"]}))
table[,"empty4"] <- NA
table[,"cohensd.ci"] <- paste(round(as.numeric(d_dx),digits=3), "(", as.character(sapply(ci_d, function(p){ paste(round(p, digits=2), collapse=",", sep="")})), ")",sep="")
table[,"BF10"] <- as.numeric(BF_BIC)
table[,"AIC"] <- as.numeric(aic_main)
table[,"BIC"] <- as.numeric(bic_main)
	

## add column to indicate whether region survives mcc
mcc <- c(table[grep("thick", table$region),][which(p.adjust(table[grep("thick|Thickness", table$region),]$p_diag, method="fdr") < 0.05),"region"],
         table[grep("surf", table$region),][which(p.adjust(table[grep("surf|SurfArea", table$region),]$p_diag, method="fdr") < 0.05),"region"],
         table[-grep("surf|thick", table$region),][which(p.adjust(table[-grep("surf|thick|SurfArea|Thickness", table$region),]$p_diag, method="fdr") < 0.05),"region"])

cond <- paste("table$region=='",mcc, "'", collapse="|", sep="")
w <- eval(parse(text=paste("which(",cond,")", sep="")))
table$mcc <- "no"
table[w,"mcc"] <- "yes"

## save table
savename <- 'Model1_results.xls'
write.table(table, file=paste(working_dir,savename, sep=""), row.names = F, col.names=T, sep="\t")


##========================================= get significant region for plot  =================================================
data.models = lapply(lme_models, function(x){ summary(x)$tTable })
p.vals = lapply(data.models, function(x){ x["Dependentanydrug1","p-value"]})
                
## apply mcc per measure
sig <- c(which(p.adjust(p.vals[grep("thick|Thickness",names(p.vals))], method="fdr") < 0.05),
         which(p.adjust(p.vals[grep("surf|SurfArea",names(p.vals))], method="fdr") < 0.05),
         which(p.adjust(p.vals[grep("surf|thick|SurfArea|Thickness",names(p.vals), invert=TRUE)], method="fdr") < 0.05) )

print(">> the following AIs showed mcc significant association with diagnosis:")
print(names(sig))



