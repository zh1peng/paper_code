# Figure 2 and its source tables.
source('lib/common.R')
run_cross_sectional <- function(input, output_dir) {
  estimates <- tests <- components <- list()
  for (co in cohorts) {
    b <- input$cross_sectional[[co]]
    b$baseline_wave <- factor(b$wave)
    mu <- mean(b$depression_score)
    sigma <- sd(b$depression_score)
    b$depression_z <- (b$depression_score - mu) / sigma
    if (co == 'NHANES') {
      # Build the design before restricting to the eligible complete-case domain.
      z <- input$nhanes_design
      z$baseline_wave <- factor(z$wave)
      z$depression_z <- (z$depression_score - mu) / sigma
      design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                          weights = ~weight_fasting, nest = TRUE, data = z)
      design <- subset(design, analysis_eligible)
      stopifnot(setequal(design$variables$id, b$id))
    }
    for (term in c('depression_binary', 'depression_z')) {
      form <- reformulate(c(term, covars, 'baseline_wave'), 'cmm')
      fit <- if (co == 'NHANES') {
        svyglm(form, design = design, family = quasibinomial())
      } else glm(form, data = b, family = binomial())
      label <- if (term == 'depression_binary') 'Primary binary' else 'Primary per SD'
      estimates[[paste(co, term)]] <- effect_row(
        fit, term, co, label, nrow(b), sum(b$cmm),
        df = if (co == 'NHANES') fit$df.residual else Inf)
      if (term == 'depression_binary') {
        if (co == 'NHANES') {
          p <- regTermTest(fit, ~depression_binary, method = 'Wald')$p
          method <- 'Survey Wald F test'
        } else {
          base <- glm(reformulate(c(covars, 'baseline_wave'), 'cmm'),
                      data = b, family = binomial())
          p <- anova(base, fit, test = 'LRT')[2, 'Pr(>Chi)']
          method <- 'Nested logistic likelihood-ratio test'
        }
        tests[[co]] <- data.frame(cohort = co, method, p = as.numeric(p))
      }
    }

    # Figure 2B appears in the main manuscript: retain its component associations.
    for (outcome in c('diabetes','heart_disease','stroke')) {
      keep <- !is.na(b[[outcome]])
      form <- reformulate(c('depression_z',covars,'baseline_wave'),outcome)
      fit <- if (co=='NHANES') {
        subdesign <- subset(design,id %in% b$id[keep])
        svyglm(form,design=subdesign,family=quasibinomial())
      } else glm(form,data=b[keep, ],family=binomial())
      components[[paste(co,outcome)]] <- effect_row(fit,'depression_z',co,outcome,
        sum(keep),sum(b[[outcome]][keep]),df=if(co=='NHANES')fit$df.residual else Inf)
    }
  }
  write_result(bind_rows(components),'component_estimates',output_dir)
  rows <- bind_rows(estimates)
  write_result(rows, 'cross_sectional_estimates', output_dir)
  write_result(bind_rows(tests), 'added_symptom_coefficient_tests', output_dir)
  pooled <- lapply(unique(rows$analysis), function(a) {
    pool_effects(rows[rows$analysis == a, ], a)
  })
  write_result(bind_rows(pooled), 'cross_sectional_meta', output_dir)
}

run_cross_sectional(load_input(),tables)
finish('2_crosssectional')
