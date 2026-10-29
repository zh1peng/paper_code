# Derive analysis visits, intervals and sample-flow counts from harmonized person-waves.
# This is analysis-sample preparation, not source-specific raw-data harmonization.
source('lib/common.R')
p <- as.data.frame(readRDS(master_file))
required <- c('cohort','id','wave','interviewed','age_eligible','depression_score',
              'depression_binary','cmm','diabetes','heart_disease','stroke',covars,
              'hypertension','time_calendar','base_count','death_time','known_dead',
              'date_conflict','CRP_raw','WBC_raw','weight_fasting','psu','stratum')
stopifnot(all(required %in% names(p)), !anyDuplicated(p[c('cohort','id','wave')]),
          all(cohorts %in% p$cohort), !anyNA(p$id), !anyNA(p$wave))
# Explicit factor levels preserve model contrasts across machines and input formats.
levels_by_field <- list(gender=c('Male','Female'),
  education_group=c('Less than High School','High School','More than High School'),
  marital_status=c('Other','Married'), smoking=c('No','Yes'), drinking=c('No','Yes'))
for (v in names(levels_by_field)) {
  stopifnot(all(is.na(p[[v]]) | as.character(p[[v]]) %in% levels_by_field[[v]]))
  p[[v]] <- factor(as.character(p[[v]]), levels=levels_by_field[[v]])
}
stopifnot(all(is.na(p$cmm) | p$cmm %in% 0:1),
          all(is.na(p$depression_binary) | p$depression_binary %in% 0:1))
p$SDMVPSU <- p$SDMVSTRA <- NA_real_
nh <- p$cohort=='NHANES'
p$SDMVPSU[nh] <- as.numeric(p$psu[nh])
p$SDMVSTRA[nh] <- as.numeric(p$stratum[nh])
first_visit <- function(d) d %>% arrange(id,wave) %>% group_by(id) %>% slice(1) %>% ungroup()
input <- list(cross_sectional=list(), longitudinal=list(), longitudinal_baseline=list(),
              endpoints=list(), biomarker=list())
flow <- cross_flow <- list()
add_flow <- function(panel, co, stages, events) {
  n <- vapply(stages, function(d) n_distinct(d$id), integer(1))
  stopifnot(all(diff(n) <= 0))
  flow[[paste(panel,co)]] <<- data.frame(panel,cohort=co,stage=seq_along(n),
    criterion=names(stages),n,excluded_from_previous=c(NA,-diff(n)),
    final_cases_or_events=c(rep(NA,length(n)-1),events))
}
for (co in cohorts) {
  cat('Preparing', co, '\n')
  panel <- p[p$cohort == co, ]
  d <- panel %>% filter(interviewed,age_eligible)
  if (co == 'NHANES') d <- d %>% filter(!is.na(weight_fasting),weight_fasting>0)
  s <- d %>% filter(!is.na(depression_score))
  candidate <- first_visit(s %>% filter(!is.na(cmm)))
  b <- candidate[complete.cases(candidate[,covars]), ]
  b$baseline_wave <- factor(b$wave)
  b$depression_z <- as.numeric(scale(b$depression_score))
  input$cross_sectional[[co]] <- b
  cross_flow[[co]] <- data.frame(cohort=co,complete=nrow(b),events=sum(b$cmm))
  add_flow('A',co,list(age_eligible=d,symptoms_available=s,CMM_available=candidate,
                       complete_covariates=b),sum(b$cmm))
  if (co == 'NHANES') {
    # Include out-of-domain records when defining strata and clusters.
    z <- panel %>% filter(!is.na(weight_fasting),weight_fasting>0)
    z$analysis_eligible <- z$id %in% b$id
    input$nhanes_design <- z
  } else {
    b0 <- first_visit(d %>% filter(!is.na(depression_score),!is.na(time_calendar),
                                   !is.na(base_count)))
    prior <- d %>% filter(cmm==1) %>% group_by(id) %>%
      summarise(first_cmm_wave=min(wave),.groups='drop')
    candidate <- left_join(b0,prior,by='id') %>%
      filter(base_count<2,is.na(first_cmm_wave)|first_cmm_wave>wave)
    b <- candidate[complete.cases(candidate[,covars]), ]
    # Standardize before excluding participants without a later evaluable visit.
    b$depression_z <- as.numeric(scale(b$depression_score))
    b$baseline_wave <- factor(b$wave)
    breaks <- unique(quantile(b$depression_score,c(0,.25,.5,.75,1)))
    b$score_group <- cut(b$depression_score,breaks=breaks,include.lowest=TRUE)
    b$start_calendar <- b$time_calendar
    # The conflict flag is participant-wide; discard the baseline row's version.
    b$date_conflict <- NULL
    records <- split(panel,panel$id)
    intervals <- endpoints <- list()
    admin <- max(panel$time_calendar,na.rm=TRUE)
    for (i in seq_len(nrow(b))) {
      bi <- b[i, ]; history <- records[[bi$id]]
      z <- history[order(history$time_calendar,history$wave), ]
      z <- z[!is.na(z$time_calendar)&!is.na(z$cmm)&
               z$time_calendar>bi$start_calendar&z$wave>bi$wave, ]
      conflict <- any(history$date_conflict,na.rm=TRUE)
      if (conflict) z <- z[FALSE, ]
      if (nrow(z)>0 && any(z$cmm==1)) z <- z[seq_len(which(z$cmm==1)[1]), ]
      if (nrow(z)>0) {
        stop <- z$time_calendar-bi$start_calendar
        start <- c(0,head(stop,-1))
        intervals[[length(intervals)+1]] <- data.frame(id=bi$id,tstart=start,
          tstop=stop,dt=stop-start,event=as.integer(z$cmm==1),end_wave=z$wave)
      }
      last <- if (nrow(z)>0) tail(z$time_calendar,1) else bi$start_calendar
      event <- as.integer(nrow(z)>0 && tail(z$cmm,1)==1)
      death <- bi$death_time
      death_event <- !conflict && event==0 && is.finite(death) && death>last && death<=admin
      status <- if (event==1) 1L else if (death_event) 2L else 0L
      end <- if (death_event) death else last
      endpoints[[i]] <- data.frame(id=bi$id,time=end-bi$start_calendar,
        observed_time=last-bi$start_calendar,event,status,date_conflict=conflict)
    }
    it <- left_join(bind_rows(intervals),b,by='id')
    it$time_mid <- (it$tstart+it$tstop)/2
    stopifnot(all(it$dt>0),all(tapply(it$event,it$id,sum)<=1))
    input$longitudinal[[co]] <- it
    input$longitudinal_baseline[[co]] <- b
    input$endpoints[[co]] <- left_join(bind_rows(endpoints),b,by='id')
    add_flow('B',co,list(age_eligible=d,score_date_diseases=b0,CMM_free=candidate,
      complete_covariates=b,valid_followup=b[b$id %in% it$id, ]),sum(it$event))
  }
  if (co != 'KLOSA') {
    blood <- d %>% filter(!is.na(CRP_raw),!is.na(WBC_raw),CRP_raw>=0,WBC_raw>0)
    candidate <- first_visit(blood %>% filter(!is.na(depression_score),!is.na(cmm)))
    b <- candidate[complete.cases(candidate[,covars]), ]
    input$biomarker[[co]] <- b
    add_flow('C',co,list(age_eligible_source=d,CRP_WBC_available=blood,
      symptoms_CMM_available=candidate,complete_covariates=b),sum(b$cmm))
  }
}
# This cache contains participant records and is LOCAL ONLY; never distribute it.
saveRDS(input,file.path(cache,'analysis_samples.rds'))
write_result(bind_rows(flow),'sample_flow_stages',tables)
write_result(bind_rows(cross_flow),'cross_sectional_flow',tables)
finish('0_prepare')
