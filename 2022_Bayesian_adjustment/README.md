
### Prepare effect sizes and sampling variances for each study

* 0_analysis_functions.R

utility functions used for effect size and sampling variance estimation

* 0_load_data.R

load freesurfer data and do necessary clean up

* 0_combat_for_ENIGMA_sMRI.R

combat function

* 1_prepare_combat_data.R

do combat for the loaded data

* 2_one_site_es_sampling_vacc_function.R

HPC function to calculate the effect sizes and sampling variances of each region for each study

### Create Hierachical models and report results
* 3_Hierachical_Bayesian_model.R

Create Hierachical models and do Jags sampling

* 4_brain_plot_and_ridge_plot_for_overarching_estimates.R

Make brain plots using ggseg and ridge plot for the overarching effect sizes

* 5_report_overarching_estimates.R

Report overarching effect sizes for each region

* 6_ridge_plot_for_each_study.R

Plot the differences between observed and posterior effect sizes for each study
