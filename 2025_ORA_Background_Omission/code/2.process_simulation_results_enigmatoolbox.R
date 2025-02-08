library(dplyr)
library(BrainEnrich)
# remotes::install_github('zh1peng/BrainEnrich')

data_path = '/media/NAS/Projects/ORA/data'
res_path='/media/NAS/Projects/ORA/sim_res_AHBA_protein_BP'

anno2test <- "GO_BP"
protein_genes <- readRDS(file.path(data_path, "protein_genes.rds")) # use use annotation genes from MF
annoData= readRDS(file.path(data_path, sprintf("%s.rds", anno2test)))
geneSetList <- get_geneSetList(annoData)
selected.gs <- filter_geneSetList(bg_genes =protein_genes, geneSetList = geneSetList, minGSSize = 20, maxGSSize = 200)
pathwaysize <- sapply(selected.gs, length)
Anno.df=Anno2Table(annoData)
pathID2Name <- setNames(Anno.df$pathName, Anno.df$pathID)

## save ahba genes for subsequent use
# for (rdonor in c("0.2", "0.4", "0.6", "0.8", "-1")) {
#   ahba_df <- read_AHBA_data(folder_path = file.path(data_path, "AHBA_ENIGMA"), rdonor = rdonor, atlas = "desikan")
#   ahba_genes <- colnames(ahba_df)
#   saveRDS(ahba_genes, sprintf('%s/AHBA_genes_rdonor%s.rds', res_path, rdonor))
# }

n_reps <- 10000

for (rdonor in c("0.2", "0.4", "0.6", "0.8", "-1")) {

ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))

for (k_perc in c(10, 20, 30)) {
  k <- round(k_perc * length(ahba_genes) / 100)
  resfile <- sprintf('%s/Simulation_data_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)

  results_df <- read.csv(resfile)
  results_df=results_df%>%
    group_by(Rep) %>%
    mutate(
      pvalsfdr.protein = p.adjust(pvals.protein, method = "fdr"),
      pvalsfdr.ahba = p.adjust(pvals.ahba, method = "fdr"),
      ifsig.protein = pvals.protein < 0.05,
      ifsig.ahba = pvals.ahba < 0.05,
      ifsigfdr.protein = pvalsfdr.protein < 0.05,
      ifsigfdr.ahba = pvalsfdr.ahba < 0.05
    ) %>%
    ungroup()

    # Calculate significant probabilities grouped by pathway
  res.psig <- results_df %>%
    group_by(Pathway) %>%
    summarise(
      psig.protein = mean(ifsig.protein, na.rm = TRUE),
      psig.ahba = mean(ifsig.ahba, na.rm = TRUE),
      psigfdr.protein = mean(ifsigfdr.protein, na.rm = TRUE),
      psigfdr.ahba = mean(ifsigfdr.ahba, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(inflation_log = log((psig.protein + 0.001)/(psig.ahba + 0.001)),  # Use log(x + 1) to avoid log(0)
    inflation_logfdr = log((psigfdr.protein + 0.001)/(psigfdr.ahba + 0.001)),
    inflation_log1 = log((psig.protein + 1)/(psig.ahba + 1)),  # Use log(x + 1) to avoid log(0)
    inflation_logfdr1 = log((psigfdr.protein + 1)/(psigfdr.ahba + 1)),  # Use log(x + 1) to avoid log(0)
    inflation_log2 = log((psig.protein + 0.1)/(psig.ahba + 0.1)),  # Use log(x + 1) to avoid log(0)
    inflation_logfdr2 = log((psigfdr.protein + 0.1)/(psigfdr.ahba + 0.1)),  # Use log(x + 1) to avoid log(0)
    PathwaySize = pathwaysize[Pathway],
    PathwayName = pathID2Name[Pathway])%>%
    arrange(desc(inflation_log))

  write.csv(res.psig, sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test), row.names = FALSE)
  }
}


