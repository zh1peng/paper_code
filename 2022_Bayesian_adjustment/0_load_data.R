#=============================
## ENIGMA ADDICTION-Bayesian Project
## Loading data
#=============================

## load packages
library(plyr)
library(dplyr) # load dplyr after plyr!!!!

print("------------------------------------")
print("Start data loading")
##----------------------------------------

## load demo data
LRdata = read.delim2(paste(working_dir,"T_final_qced_30days.txt", sep=""),sep = ',')

## check column classes
table(apply(LRdata, 2, class)) ## all columns are of class character

## set LR data to numeric
LRdata2 <- LRdata
for(col in colnames(LRdata2)[grep("L_|R_|Age|Alclast30days|Niclast30days|ICV", colnames(LRdata))] ){
  LRdata2[,col] <- as.character(LRdata[,col])
  LRdata2[,col] <- as.numeric(LRdata2[,col])
}
TRUE %in% apply(LRdata==LRdata2, 2, function(x){FALSE %in% x}) ## checked, values don't change after setting to numeric
LRdata <- LRdata2
rm(col, LRdata2)


## note: reset levels for factor
LRdata$Site <- as.character(LRdata$Site)
LRdata$Site <- as.factor(LRdata$Site)
#LRdata$Age <- as.numeric(as.character(LRdata$Age)) #caution here!!!
#=========================================
## Define reference groups
LRdata$Dependentanydrug <- as.factor(as.character(LRdata$Dependentanydrug))
LRdata$Dependentanydrug <- relevel(LRdata$Dependentanydrug, "0")

LRdata$DependentonPrimaryDrug <- as.factor(as.character(LRdata$DependentonPrimaryDrug))
LRdata$DependentonPrimaryDrug <- relevel(LRdata$DependentonPrimaryDrug,"0")

LRdata$PrimaryDrug <- as.factor(as.character(LRdata$PrimaryDrug))
LRdata$PrimaryDrug <- relevel(LRdata$PrimaryDrug,"0")

LRdata$Sex <- as.factor(LRdata$Sex)
LRdata$Sex <- relevel(LRdata$Sex, "1")



##---------------------------------
print("remove NULL values from LR data")
##---------------------------------
## remove all the NULL values from the LR data columns
## loop through all regions
for(i in colnames(LRdata)[grep("L_|R_", colnames(LRdata))] ){
  if(length(which(LRdata[,i]==0)) !=0){
    print(paste0("remove ",length(which(LRdata[,i]==0))," NULLs for ",i))
    LRdata[which(LRdata[,i]==0),i] <- NA
    LRdata[,i]=as.numeric(LRdata[,i])
  }
}



if(remove_dataset=="yes"){
  #-------------------------------------------------------------
  print("remove datasets with case-only or control-only data")
  #------------------------------------------------------------
  #library(dplyr)
  tmp = LRdata %>% group_by(DependentonPrimaryDrug, Site) %>% count(Dependentanydrug)
  
  if( length(names(which(table(tmp$Site)==1))) != 0) ## if 1 instead of 2, only cases or only control data is present
  {
    only = names(which(table(tmp$Site)==1))
    cond = paste("LRdata$Site=='",only,"'", sep="", collapse="|")
    w <- eval(parse(text=paste("which(", cond, ")")))
    LRdata <- LRdata[-w,]; rm(cond, w, tmp)
    
    print(" >> The following datasets were removed:"); print(only)
    
  }
}


covs <- c("PI","Subject","Dependentanydrug","DependentonPrimaryDrug", "PrimaryDrug","Age", "Sex", "Site","Half","Alclast30days","Niclast30days","ICV")
covs=c(covs, colnames(LRdata)[grep("L_|R_", colnames(LRdata))] )

LR2use=LRdata[,covs]

print("finished data loading")
print("------------------------------------")

####
##
#
