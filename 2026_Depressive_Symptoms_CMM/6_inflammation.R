# Figure 6: joint-category prevalence ratios and main additive contrasts.
source('lib/common.R')
run_biomarker <- function(input, output_dir) {
  groups <- additive <- list()
  for (co in c('NHANES', 'CHARLS', 'HRS')) {
    d <- input$biomarker[[co]]
    stopifnot(!anyDuplicated(d$id), all(d$CRP_raw >= 0), all(d$WBC_raw > 0))
    d$baseline_wave <- factor(d$wave)
    score <- as.numeric(scale(d$CRP_raw)) + as.numeric(scale(d$WBC_raw))
    d$high <- as.integer(score >= median(score))
    terms <- c('depression_binary*high', covars)
    if (nlevels(d$baseline_wave) > 1) terms <- c(terms, 'baseline_wave')
    form <- reformulate(terms, 'cmm')
    if (co == 'NHANES') {
      design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                          weights = ~weight_fasting, nest = TRUE, data = d)
      fit <- svyglm(form, design = design, family = quasipoisson(link = 'log'))
      V <- vcov(fit)
      df <- fit$df.residual
    } else {
      fit <- glm(form, data = d, family = poisson(link = 'log'))
      V <- vcovHC(fit, type = 'HC0')
      df <- Inf
    }
    terms <- c('depression_binary', 'high', 'depression_binary:high')
    beta <- coef(fit)[terms]
    V <- V[terms, terms]
    q <- qt(.975, df)
    for (a in 0:1) for (b in 0:1) {
      contrast <- c(a, b, a*b)
      eta <- sum(contrast * beta)
      se <- sqrt(as.numeric(t(contrast) %*% V %*% contrast))
      keep <- d$depression_binary == a & d$high == b
      groups[[paste(co, a, b)]] <- data.frame(
        cohort = co, analysis = 'Primary', depressive_symptoms = a, high_inflammation = b,
        n = sum(keep), cases = sum(d$cmm[keep]), PR = exp(eta),
        lcl = exp(eta - q*se), ucl = exp(eta + q*se))
    }
    pr <- c(exp(beta[1]), exp(beta[2]), exp(sum(beta)))
    value <- unname(pr[3] - pr[1] - pr[2] + 1)
    gradient <- c(pr[3] - pr[1], pr[3] - pr[2], pr[3])
    se <- sqrt(as.numeric(t(gradient) %*% V %*% gradient))
    additive[[co]] <- data.frame(cohort = co, contrast = value,
                                 lcl = value - q*se, ucl = value + q*se)
  }
  write_result(bind_rows(groups), 'biomarker_groups', output_dir)
  write_result(bind_rows(additive), 'biomarker_additive', output_dir)
}

run_biomarker(load_input(),tables)
finish('6_inflammation')
