# Figure 5: descriptive Aalen-Johansen CMM incidence, with recorded death competing.
# This is not Kaplan-Meier risk and does not add a supplementary Cox model.
source('lib/common.R')
input <- load_input()
curves <- risk <- list()
for (co in cohorts[-1]) {
  e <- input$endpoints[[co]] %>%
    filter(time>0,!date_conflict,!(known_dead & is.na(death_time) & event==0))
  for (g in levels(e$score_group)) {
    z <- e[e$score_group==g, ]
    z$state <- factor(z$status,levels=0:2,labels=c('censored','CMM','death'))
    f <- survfit(Surv(time,state)~1,data=z)
    curves[[paste(co,g)]] <- data.frame(cohort=co,group=g,time=c(0,f$time),
      CIF=c(0,f$pstate[,'CMM']),n_risk=c(nrow(z),f$n.risk[,1]),initial_n=nrow(z))
    times <- seq(0,floor(max(e$time)/2)*2,2)
    risk[[paste(co,g)]] <- data.frame(cohort=co,group=g,time=times,
      n=sapply(times,function(t)sum(z$time>=t)))
  }
}
write_result(bind_rows(curves),'descriptive_cif',tables)
write_result(bind_rows(risk),'cif_at_risk',tables)
finish('5_survival')
