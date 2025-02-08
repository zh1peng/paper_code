library(ggplot2)
library(patchwork)

data_path='E:/xhmhc/ORA misuse/data'

bg_type <- "protein"  # or annotation
anno2test <- "GO_BP"


plots <- list()
# Create a list to store plots

# Iterate through donors and percentage values
  for (k_perc in c(10, 20, 30)) {
for (anno2test in c('GO_BP','GO_MF')){
for (bg_type in c('annotation','protein')){
  
 res_path = sprintf('E:/xhmhc/ORA misuse/sim_res_Aurina_%s_%s', bg_type, substr(anno2test, 4, 5))

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
rdonor = 'no'


  ahba_genes <- readRDS(sprintf('%s/AHBA_genes_Aurina2019.rds', data_path))
    k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))
    # Determine if the plot is in the first column or last row
    # show_x_axis <- (k_perc == 30)  # Show x-axis only for the last row (k_perc == 30)
    show_x_axis <- (k_perc == 30)
    show_y_axis <- (anno2test == "GO_BP" & bg_type =="annotation")  # Show y-axis only for the first column (rdonor == "0.2")
    show_title <- (k_perc == 10)  # Show title only for the first column (rdonor == "0.2")
    # Generate the scatter plot
    p.psig <- ggplot(res.psig, aes(y = psig.protein, x = psig.ahba)) +
      geom_point(aes(color = ratio), size = 1.5, alpha = 0.8) +
      scale_color_viridis_c() +
      annotate("text", x = 0.5, y = 1, label = sprintf("k=%d (%d%%)", k, k_perc), 
               hjust = 0.5, size = 4, color = "black") +
      labs(
        y = if (show_y_axis) sprintf("%s (n=%d)", anno_text, length(protein_genes)) else NULL,
        x = if (show_x_axis) sprintf("AHBA (n=%d)", length(ahba_genes)) else NULL,
        color = "Ratio",
        title = if (show_title) sprintf("%s %s",  substr(anno2test, 4, 5), anno_text) else NULL
      ) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = c(0.80, 0.36),
        legend.title = element_text(size = 12),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.4, "cm"),

        axis.title.x = if (show_x_axis) element_text() else element_blank(),
        axis.title.y = if (show_y_axis) element_text() else element_blank(),
        axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
        axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
        axis.text.x = if (show_x_axis) element_text(size = 12) else element_blank(),
        axis.text.y = if (show_y_axis) element_text(size = 12) else element_blank(),
        plot.title = if (show_title) element_text(hjust = 0.5,size=12) else element_blank() 
      ) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#3f3f3f") +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "#3f3f3f") +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "#3f3f3f") +
      coord_fixed(1, xlim = c(0, 1), ylim = c(0, 1))

    # Add the plot to the list
    plots[[paste(anno2test, bg_type, k_perc)]] <- p.psig
  }
  }
}
 final_p=wrap_plots(plots, nrow=3, ncol = 4) + 
 plot_annotation(title = sprintf("Background: AHBA vs. Others"),
  theme = theme(plot.title = element_text(size = 14, hjust = 0.5,face = "bold")))

ggsave(sprintf('%s/final_full_scatter.png', res_path),final_p, width = 11, height = 9, dpi = 300)

