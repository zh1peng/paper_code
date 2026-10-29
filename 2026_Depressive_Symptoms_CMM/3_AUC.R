# Figure 3: nested logistic models, paired refitted bootstrap, optimism and calibration.
# This reproduces the manuscript's discrimination analysis; no extra ML algorithms.
source('lib/common.R')
weighted_auc <- function(y, p, w) {
  order <- order(p)
  y <- y[order]; p <- p[order]; w <- w[order]
  group <- cumsum(c(TRUE, diff(p) != 0))
  pos <- as.numeric(rowsum(w*y, group))
  neg <- as.numeric(rowsum(w*(1-y), group))
  sum(pos * (cumsum(neg) - neg/2)) / (sum(pos) * sum(neg))
}

run_discrimination <- function(baselines, output_dir, B = 5000) {
  results <- list()
  for (co in cohorts) {
    set.seed(21506 + match(co, cohorts))
    d <- baselines[[co]]
    d$baseline_wave <- factor(d$wave)
    X0 <- model.matrix(reformulate(c(covars, 'baseline_wave')), d)
    X1 <- cbind(X0, depression_binary = d$depression_binary)
    y <- d$cmm
    w <- if (co == 'NHANES') d$weight_fasting else rep(1, nrow(d))
    w <- w/mean(w)
    f0 <- glm.fit(X0, y, weights = w, family = quasibinomial())
    f1 <- glm.fit(X1, y, weights = w, family = quasibinomial())
    a0 <- weighted_auc(y, plogis(X0 %*% coef(f0)), w)
    a1 <- weighted_auc(y, plogis(X1 %*% coef(f1)), w)
    if (co == 'NHANES') {
      design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                          weights = ~weight_fasting, nest = TRUE, data = d)
      replicate <- bootweights(design$strata[,1], design$cluster[,1],
                                replicates = B, compress = TRUE)
      rw <- as.matrix(replicate$repweights)/design$prob
    }
    boot <- matrix(NA_real_, B, 5,
                   dimnames = list(NULL, c('delta', 'optimism0', 'optimism1',
                                           'test_slope', 'test_intercept')))
    for (i in seq_len(B)) {
      wb <- if (co == 'NHANES') rw[,i] else {
        tabulate(sample.int(nrow(d), nrow(d), replace = TRUE), nbins = nrow(d))
      }
      wb <- wb/mean(wb)
      b0 <- glm.fit(X0, y, weights = wb, family = quasibinomial())
      b1 <- glm.fit(X1, y, weights = wb, family = quasibinomial())
      stopifnot(b0$converged, b1$converged,
                all(is.finite(coef(b0))), all(is.finite(coef(b1))))
      eta <- as.numeric(X1 %*% coef(b1))
      p0 <- plogis(X0 %*% coef(b0)); p1 <- plogis(eta)
      train0 <- weighted_auc(y, p0, wb); train1 <- weighted_auc(y, p1, wb)
      test0 <- weighted_auc(y, p0, w); test1 <- weighted_auc(y, p1, w)
      calibration <- glm.fit(cbind(1, eta), y, weights = w, family = quasibinomial())
      intercept <- uniroot(function(a) sum(w*(y-plogis(a+eta))), c(-20,20))$root
      boot[i,] <- c(train1-train0, train0-test0, train1-test1,
                    coef(calibration)[2], intercept)
      if (i %% 250 == 0) cat(co, i, '/', B, 'bootstrap replicates\n')
    }
    stopifnot(all(is.finite(boot)))
    ci <- quantile(boot[,'delta'], c(.025,.975))
    results[[co]] <- data.frame(
      cohort = co, n = nrow(d), AUC_core = a0, AUC_added = a1, delta = a1-a0,
      delta_lcl = unname(ci[1]), delta_ucl = unname(ci[2]),
      corrected_AUC_core = a0-mean(boot[,'optimism0']),
      corrected_AUC_added = a1-mean(boot[,'optimism1']),
      corrected_calibration_slope = mean(boot[,'test_slope']),
      corrected_calibration_intercept = mean(boot[,'test_intercept']), replicates = B)
    saveRDS(boot, file.path(cache, paste0(co, '_auc_bootstrap.rds')))
  }
  # Percentile intervals describe apparent paired AUC differences, not a null test.
  write_result(bind_rows(results), 'discrimination_bootstrap', output_dir)
}


B <- if (length(args)>=3) as.integer(args[3]) else 5000L
stopifnot(!is.na(B),B>=2)
run_discrimination(load_input()$cross_sectional,tables,B=B)
finish('3_AUC')
