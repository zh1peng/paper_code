packages <- c('dplyr','survey','metafor','sandwich','survival')
missing <- setdiff(packages,rownames(installed.packages()))
if (length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
