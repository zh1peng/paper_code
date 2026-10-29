# Core model definitions used by the revised CMM manuscript.
suppressPackageStartupMessages({
  library(dplyr)
  library(survey)
  library(metafor)
  library(sandwich)
  library(survival)
})
options(survey.lonely.psu = 'adjust', survey.adjust.domain.lonely = TRUE)
cohorts <- c('NHANES', 'CHARLS', 'HRS', 'KLOSA')
covars <- c('age', 'gender', 'education_group', 'marital_status',
            'smoking', 'drinking', 'bmi')

effect_row <- function(fit, term, cohort, analysis, n, events,
                       V = vcov(fit), df = Inf) {
  b <- unname(coef(fit)[term])
  se <- unname(sqrt(diag(V))[term])
  q <- qt(.975, df)
  data.frame(cohort, analysis, n, events, estimate = b, se,
             p = 2 * pt(-abs(b/se), df), ratio = exp(b), lcl = exp(b - q * se), ucl = exp(b + q * se))
}

pool_effects <- function(rows, analysis) {
  f <- rma.uni(yi = rows$estimate, sei = rows$se, method = 'REML', test = 'knha')
  b <- as.numeric(f$b)
  s <- sqrt(f$tau2 + f$se^2)
  k <- nrow(rows)
  data.frame(analysis, k, ratio = exp(b), lcl = exp(f$ci.lb), ucl = exp(f$ci.ub),
             p = f$pval, tau2 = f$tau2, I2 = f$I2, Q = f$QE, Q_p = f$QEp,
             PI_k1_lower = exp(b - qt(.975, k - 1) * s),
             PI_k1_upper = exp(b + qt(.975, k - 1) * s),
             PI_k2_lower = exp(b - qt(.975, k - 2) * s),
             PI_k2_upper = exp(b + qt(.975, k - 2) * s))
}

write_result <- function(x, name, output_dir) {
  write.csv(x, file.path(output_dir, paste0(name, '.csv')), row.names = FALSE, na = '')
}


# Every numbered script uses the same command-line contract.
# Run from this package directory: Rscript SCRIPT.R master_person_wave.rds output
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Supply preprocessed master RDS and output directory')
master_file <- args[1]
output_root <- args[2]
tables <- file.path(output_root, 'tables')
cache <- file.path(output_root, 'cache')
for (path in c(tables, cache, file.path(output_root, 'figures'))) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
load_input <- function() readRDS(file.path(cache, 'analysis_samples.rds'))
finish <- function(stage) {
  capture.output(sessionInfo(), file = file.path(output_root, paste0('session_', stage, '.txt')))
}
