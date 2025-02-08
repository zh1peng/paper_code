library(dplyr)
library(ggplot2)
library(patchwork)

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


plist <- list()

# make bar plots
for (k_perc in c(10,20,30)) {
for (rdonor in c("0.6","0.4","0.2", "-1")) {
ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))

  
top_pathways <- res.psig %>%
  arrange(desc(psig.protein)) %>%
  head(10) %>% # Select the top 10 pathways
  mutate(Pathway = forcats::fct_reorder(Pathway, psig.protein))
  top_pathways$PathwayName <- stringr:: str_wrap(top_pathways$PathwayName, width = 30)  # Adjust width as needed

show_x_axis <- TRUE  # Show x-axis only for the last row (k_perc == 30)
show_x_label <- (k_perc == 30)  # Show x-axis label only for the last row (k_perc == 30)
 # Create the bar plot
p <- ggplot(top_pathways, aes(x = psig.protein, y = Pathway, fill = ratio)) +
  geom_bar(stat = "identity", width = 0.8) +  # Bar plot
  scale_fill_viridis_c(option = "viridis", limits=range(res.psig$ratio)) +  # Color scale for psig.protein
  geom_text(aes(label = PathwayName),  # Add pathway names above bars
            hjust = 1.05, vjust = 0.5, size = 2.5, color = "black",lineheight=0.8) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +  # Wrap long pathway names
  theme_bw(base_size = 14) +  # Clean layout
  theme(
    axis.title.y = element_blank(),  # No y-axis label
    axis.title.x = if (show_x_label) element_text(size = 12) else element_blank(),  # X-axis text size
    axis.text.y = element_blank(),  # No y-axis text
    axis.text.x = if (show_x_axis) element_text(size = 10) else element_blank(),  # X-axis text size
    axis.ticks.x = if (show_x_axis) element_line() else element_blank(),  # X-axis ticks
    axis.ticks.y = element_blank(),  # No y-axis ticks
    plot.title = element_text(size = 12, hjust = 0.5),  # Title
    legend.position = "none"  # No legend
  ) +
  labs(
    title = sprintf("AHBA genes (n=%d); k=%d (%d%%)", length(ahba_genes), k, k_perc),
    x = sprintf("Probability of significance")
  )
 plist[[paste(rdonor, k_perc)]] <- p
}
}

 final_p=wrap_plots(plist, nrow=3, ncol = 4) + 
 plot_annotation(title = sprintf("Top 10 Pathways for %s background (n=%d)", anno_text, length(protein_genes)),
  theme = theme(plot.title = element_text(size = 14, hjust = 0.5,face = "bold")))

ggsave(sprintf('%s/full_Top10_pathways.png', res_path), plot = final_p, width =15, height =10, dpi = 300)






plist <- list()

# make bar plots
for (k_perc in c(20)) {
for (rdonor in c("0.6","0.4","0.2")) {
ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))


  
top_pathways <- res.psig %>%
  arrange(desc(psig.protein)) %>%
  head(10) %>% # Select the top 10 pathways
  mutate(Pathway = forcats::fct_reorder(Pathway, psig.protein))
top_pathways$PathwayName <- stringr:: str_wrap(top_pathways$PathwayName, width = 35)  # Adjust width as needed

show_x_axis <- (k_perc == 20)  # Show x-axis only for the last row (k_perc == 30)
 # Create the bar plot
p <- ggplot(top_pathways, aes(x = psig.protein, y = Pathway, fill = ratio)) +
  geom_bar(stat = "identity", width = 0.8) +  # Bar plot
  scale_fill_viridis_c(option = "viridis", limits=range(res.psig$ratio)) +  # Color scale for psig.protein
  geom_text(aes(label = PathwayName),  # Add pathway names above bars
            hjust = 1.05, vjust = 0.5, size = 2.5, color = "black",lineheight = 0.8) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +  # Wrap long pathway names
  theme_bw(base_size = 14) +  # Clean layout
  theme(
    axis.title.y = element_blank(),  # No y-axis label
    axis.title.x = if (show_x_axis) element_text(size = 12) else element_blank(),  # X-axis text size
    axis.text.y = element_blank(),  # No y-axis text
    axis.text.x = if (show_x_axis) element_text(size = 8) else element_blank(),  # X-axis text size
    axis.ticks.x = if (show_x_axis) element_line() else element_blank(),  # X-axis ticks
    axis.ticks.y = element_blank(),  # No y-axis ticks
    plot.title = element_text(size = 12, hjust = 0.5),  # Title
    legend.position = "none"  # No legend
  ) +
#   labs(
#     title = sprintf("AHBA (n=%d); k=%d (%d%%)", length(ahba_genes), k, k_perc),
#     x = sprintf("Probability of significance")
#   )
   labs(
    title = NULL,
    x = NULL
  )
 plist[[paste(rdonor, k_perc)]] <- p
}
}

 final_p=wrap_plots(plist, nrow=1, ncol = 3) 
#  plot_annotation(title = sprintf("Top 10 Pathways for %s background (n=%d)", anno_text, length(protein_genes)),
#   theme = theme(plot.title = element_text(size = 14, hjust = 0.5,face = "bold")))

ggsave(sprintf('%s/selected_Top10_pathways.png', res_path), plot = final_p, width =8, height =2.6, dpi = 300)
