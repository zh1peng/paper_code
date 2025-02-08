library(parallel)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(BrainEnrich)

code_path ='/media/NAS/Projects/ORA/code'
data_path = '/media/NAS/Projects/ORA/data'
res_path='/media/NAS/Projects/ORA/sim_res_AHBA_annotation_BP'
if (!dir.exists(res_path)) {
  dir.create(res_path)
}

data_path='/media/NAS/Projects/ORA/data'

source(sprintf('%s/functions.R',code_path))
sourceCpp(sprintf('%s/intersectToList.cpp',code_path))

anno2test <- "GO_BP"
annoData=readRDS(file.path(data_path, sprintf("%s.rds", anno2test)))
geneSetList <- get_geneSetList(annoData)
protein_genes <- readRDS(file.path(data_path, "bp_genes.rds"))
selected.gs <- filter_geneSetList(bg_genes =protein_genes, geneSetList = geneSetList, minGSSize = 20, maxGSSize = 200)

for (rdonor in c("0.2","0.4","0.6", "0.8", "-1")) {
ahba_genes <- readRDS(file.path(data_path, sprintf("AHBA_genes_rdonor%s.rds", rdonor)))

# Number of repetitions
n_reps <- 10000
n_cores <- detectCores() - 1 

for (k_perc in c(10, 20, 30)) {
# Number of available cores
# Create a cluster
k <- round(k_perc * length(ahba_genes) / 100)
cl <- makeCluster(n_cores)
# Export necessary variables to the cluster
clusterExport(cl, varlist = c("ahba_genes", "k", "selected.gs", "protein_genes", "fast_ORA","code_path","simulate_repetition"), envir = .GlobalEnv)
clusterEvalQ(cl, {
  Rcpp::sourceCpp(sprintf('%s/intersectToList.cpp',code_path))# Ensure the function is loaded in workers
})
results <- parLapply(cl, 1:n_reps, function(rep_id) {
  simulate_repetition(rep_id, k, selected.gs, ahba_genes,protein_genes)
})
stopCluster(cl)
results_df <- do.call(rbind, results)
savefile <- sprintf('%s/Simulation_data_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
write.csv(results_df, file = savefile, row.names = FALSE)
}
}


