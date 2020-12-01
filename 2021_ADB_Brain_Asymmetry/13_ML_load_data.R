#=============================
#Load data and calculate AI
#aidata is used for further analysis

library(dplyr)
##===========================load data ===================================================
LRdata = read.delim2(paste(working_dir,"T_final_qced_30days.txt", sep=""),sep = ',')
table(apply(LRdata, 2, class)) ## all columns are of class character

##===========================Set numeric cols==============================================
LRdata2 <- LRdata
for(col in colnames(LRdata2)[grep("L_|R_|Age|Alclast30days|Niclast30days", colnames(LRdata))] ){
  LRdata2[,col] <- as.character(LRdata[,col])
  LRdata2[,col] <- as.numeric(LRdata2[,col])
}
TRUE %in% apply(LRdata==LRdata2, 2, function(x){FALSE %in% x}) ## checked, values don't change after setting to numeric
LRdata <- LRdata2
rm(col, LRdata2)


##===========================Set factor cols==================================================
LRdata$Dependentanydrug <- as.factor(as.character(LRdata$Dependentanydrug))
LRdata$Dependentanydrug <- relevel(LRdata$Dependentanydrug, "0")

LRdata$DependentonPrimaryDrug <- as.factor(as.character(LRdata$DependentonPrimaryDrug))
LRdata$DependentonPrimaryDrug <- relevel(LRdata$DependentonPrimaryDrug,"0")

LRdata$PrimaryDrug <- as.factor(as.character(LRdata$PrimaryDrug))
LRdata$PrimaryDrug <- relevel(LRdata$PrimaryDrug,"0")

LRdata$Sex <- as.factor(LRdata$Sex)
LRdata$Sex <- relevel(LRdata$Sex, "1")

LRdata$Site <- as.character(LRdata$Site)
LRdata$Site <- as.factor(LRdata$Site)



##=============================remove Null values============================================

for(i in colnames(LRdata)[grep("L_|R_", colnames(LRdata))] ){
  if(length(which(LRdata[,i]==0)) !=0){
    print(paste0("remove ",length(which(LRdata[,i]==0))," NULLs for ",i))
    LRdata[which(LRdata[,i]==0),i] <- NA
    LRdata[,i]=as.numeric(LRdata[,i])
  }
}


##=============================calculate asymmetry index ====================================
for (i in colnames(LRdata)[grep("L_",colnames(LRdata))] )
{
  left <-  i
  right <- gsub("L_", "R_", i)
  AI <-   gsub("L_", "AI_", i)
  LRdata[,AI] <- (LRdata[,left] - LRdata[,right]) / (LRdata[,left] + LRdata[,right]) 
  rm(left, right, AI)
} 
rm(i) 


##================================remove 3std outliers =========================================
for(region in colnames(LRdata)[grep("AI_", colnames(LRdata))]){
    wins3.low=mean(LRdata[,region],na.rm = TRUE)-3*as.numeric(sd(LRdata[,region],na.rm = TRUE))
    wins3.up=mean(LRdata[,region],na.rm = TRUE)+3*as.numeric(sd(LRdata[,region],na.rm = TRUE))
    LRdata[which(LRdata[,region]>wins3.up),region] <- wins3.up
    LRdata[which(LRdata[,region]<wins3.low),region] <- wins3.low
}

##=================== remove datasets with case-only or control-only data =====================
tmp = LRdata %>% group_by(Dependentanydrug, Site) %>% count(Dependentanydrug)
if( length(names(which(table(tmp$Site)==1))) != 0) ## if 1 instead of 2, only cases or only control data is present
  {
    only = names(which(table(tmp$Site)==1))
    cond = paste("LRdata$Site=='",only,"'", sep="", collapse="|")
    w <- eval(parse(text=paste("which(", cond, ")")))
    LRdata <- LRdata[-w,]; rm(cond, w, tmp)
    
   print(" >> The following datasets were removed:"); print(only)

  }

##==============================select cols for the further analysis=================================
covs <- c("Subject", "Dependentanydrug","DependentonPrimaryDrug", "PrimaryDrug","Age", "Sex", "Site","Half","Alclast30days","Niclast30days")
aidata <- LRdata[,c(covs, colnames(LRdata)[grep("AI_", colnames(LRdata))])]






