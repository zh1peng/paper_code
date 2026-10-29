# Table 1: unweighted descriptive statistics and signed standardized differences.
# No P values. CMM component counts define the outcome, so their SMDs are omitted.
source('lib/common.R')
input <- load_input()
rows <- list()
for (co in cohorts) {
  b <- input$cross_sectional[[co]]
  b$component_count <- factor(rowSums(b[,c('diabetes','heart_disease','stroke')]),levels=0:3)
  b$elevated_symptoms <- factor(b$depression_binary,levels=0:1,labels=c('No','Yes'))
  b$hypertension <- factor(b$hypertension,levels=0:1,labels=c('No','Yes'))
  for (v in c(covars,'depression_score','elevated_symptoms','component_count','hypertension')) {
    a <- b[[v]][b$cmm==0]; z <- b[[v]][b$cmm==1]
    if (is.numeric(b[[v]])) {
      smd <- (mean(z,na.rm=TRUE)-mean(a,na.rm=TRUE))/sqrt((var(a,na.rm=TRUE)+var(z,na.rm=TRUE))/2)
      rows[[length(rows)+1]] <- data.frame(cohort=co,variable=v,level='',
        without_CMM=sprintf('%.2f (%.2f)',mean(a,na.rm=TRUE),sd(a,na.rm=TRUE)),
        with_CMM=sprintf('%.2f (%.2f)',mean(z,na.rm=TRUE),sd(z,na.rm=TRUE)),SMD=smd)
    } else for (lv in levels(b[[v]])) {
      p0 <- mean(a==lv,na.rm=TRUE); p1 <- mean(z==lv,na.rm=TRUE)
      smd <- if (v=='component_count') NA_real_ else
        (p1-p0)/sqrt((p0*(1-p0)+p1*(1-p1))/2)
      rows[[length(rows)+1]] <- data.frame(cohort=co,variable=v,level=lv,
        without_CMM=sprintf('%d (%.1f%%)',sum(a==lv,na.rm=TRUE),100*p0),
        with_CMM=sprintf('%d (%.1f%%)',sum(z==lv,na.rm=TRUE),100*p1),SMD=smd)
    }
    if (anyNA(c(a,z))) rows[[length(rows)+1]] <- data.frame(cohort=co,variable=v,
      level='Missing',without_CMM=as.character(sum(is.na(a))),
      with_CMM=as.character(sum(is.na(z))),SMD=NA_real_)
  }
}
write_result(bind_rows(rows),'table1',tables)
finish('1_summary')
