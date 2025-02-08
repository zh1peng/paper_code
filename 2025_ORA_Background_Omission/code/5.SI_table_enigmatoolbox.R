library(dplyr)
data_path='E:/xhmhc/ORA misuse/data'

bg_type <- "protein"  # or annotation
anno2test <- "GO_BP"

res_path = sprintf('E:/xhmhc/ORA misuse/sim_res_AHBA_%s_%s', bg_type, substr(anno2test, 4, 5))

if (bg_type == "protein") {
  protein_genes <- readRDS(file.path(data_path, "protein_genes.rds"))
}else if (bg_type == "annotation") {
   protein_genes <- readRDS(file.path(data_path, sprintf("%s_genes.rds", tolower(substr(anno2test, 4, 5)))))
} else {
  stop("Invalid background type. Choose 'protein' or 'annotation'.")
}
anno_text <- if (bg_type == "protein") "Protein-coding" else "Annotation"

n_reps=10000
epsilon=1/n_reps



# for (k_perc in c(10, 20)) {
#     for (rdonor in c("0.6","0.4", "0.2")) {
#   ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
#     k <- round(k_perc * length(ahba_genes) / 100)
#     res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
#     res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))
#     res2save <- res.psig %>% select(Pathway, PathwayName, PathwaySize, psig.protein, psig.ahba, ratio) %>% arrange(desc(psig.protein))
#     write.csv(res2save, sprintf('%s/clean_Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test), row.names = FALSE)
#     }
# }


library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Loop through `k_perc` and `rdonor` combinations
for (k_perc in c(10, 20)) {
  for (rdonor in c("0.6", "0.4", "0.2")) {
    ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
    k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', 
                            res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% 
                mutate(ratio = (psig.protein + epsilon) / (psig.ahba + epsilon))
    
    # Clean and arrange the data
    res2save <- res.psig %>% 
                select(Pathway, PathwayName, PathwaySize, psig.protein, psig.ahba, ratio) %>% 
                arrange(desc(psig.protein))
    
    # Define sheet name
    sheet_name <- sprintf("k%d_r%s", k_perc, rdonor)
    
    # Add the data to the workbook
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, res2save)
  }
}

# Save the workbook to a file
saveWorkbook(wb, file = sprintf('%s/Clean_Psig_Results.xlsx', res_path), overwrite = TRUE)
