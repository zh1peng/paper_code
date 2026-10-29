# Figure 4A: incident CMM in observed assessment intervals; Figure 4B: score splines.
source('lib/common.R')
input <- load_input()
effects <- curves <- followup <- list()
for (co in cohorts[-1]) {
  it <- input$longitudinal[[co]]
  b <- input$longitudinal_baseline[[co]]
  rhs <- c('depression_z',covars,'factor(base_count)','baseline_wave',
           'splines::ns(log1p(time_mid),df=3)','offset(log(dt))')
  fit <- glm(reformulate(rhs,'event'),data=it,family=binomial(link='cloglog'))
  stopifnot(fit$converged,all(is.finite(coef(fit))))
  V <- vcovCL(fit,cluster=it$id,type='HC0')
  effects[[co]] <- effect_row(fit,'depression_z',co,'Discrete-time primary',
                             n_distinct(it$id),sum(it$event),V)
  # Quantiles and reference use the full complete baseline, before follow-up exclusion.
  knots <- as.numeric(quantile(b$depression_score,c(.1,.5,.9)))
  basis <- splines::ns(it$depression_score,knots=knots[2],Boundary.knots=knots[c(1,3)])
  it$s1 <- basis[,1]; it$s2 <- basis[,2]
  sf <- glm(reformulate(c('s1','s2',setdiff(rhs,'depression_z')),'event'),
            data=it,family=binomial(link='cloglog'))
  stopifnot(sf$converged)
  sv <- vcovCL(sf,cluster=it$id,type='HC0')[c('s1','s2'),c('s1','s2')]
  grid <- seq(min(b$depression_score),max(b$depression_score),length.out=101)
  contrast <- sweep(predict(basis,grid),2,as.numeric(predict(basis,median(b$depression_score))))
  eta <- as.numeric(contrast %*% coef(sf)[c('s1','s2')])
  se <- sqrt(rowSums((contrast %*% sv)*contrast))
  curves[[co]] <- data.frame(cohort=co,score=grid,ratio=exp(eta),
    lcl=exp(eta-1.96*se),ucl=exp(eta+1.96*se),knot10=knots[1],knot50=knots[2],
    knot90=knots[3],reference=median(b$depression_score))
  ep <- input$endpoints[[co]]
  observed <- ep$observed_time[ep$observed_time>0 & !ep$date_conflict]
  followup[[co]] <- data.frame(cohort=co,n=n_distinct(it$id),events=sum(it$event),
    intervals=nrow(it),median_followup=median(observed),max_followup=max(observed))
}
rows <- bind_rows(effects)
write_result(rows,'longitudinal_estimates',tables)
write_result(pool_effects(rows,'Discrete-time primary'),'longitudinal_meta',tables)
write_result(bind_rows(curves),'longitudinal_splines',tables)
write_result(bind_rows(followup),'longitudinal_followup',tables)
finish('4_longitudinal')
