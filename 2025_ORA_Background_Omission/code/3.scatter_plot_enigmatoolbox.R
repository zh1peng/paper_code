library(ggplot2)
library(patchwork)

data_path='E:/xhmhc/ORA misuse/data'

bg_type <- "annotation"  # or annotation
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


# for (rdonor in c("0.2", "0.4", "0.6", "0.8", "-1")) {
# ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', res_path, rdonor))
# for (k_perc in c(10, 20, 30)) {
#   k <- round(k_perc * length(ahba_genes) / 100)
#   res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
#   res.psig <- read.csv(res.sigfile)
#   p.psig=ggplot(res.psig, aes(y = psig.protein, x = psig.ahba)) +
#   geom_point(aes(color = inflation_log1), size = 1.5, alpha = 0.8) +  # Pathway as color
#   scale_color_viridis_c()+
  
#    annotate("text", x = 0.5, y = 1, label = sprintf("Subsample: k=%d (%d%%)", k,k_perc), 
#            hjust = 0.5, size = 4,  color = "black")+
#   labs(
#     title = sprintf("Probability of significance (k=%d)", k),
#     #subtitle = sprintf("Number of significant genes: %d", k),
#     y = sprintf("Protein-coding genes (n=%d)" , length(protein_genes)),
#     x = sprintf("AHBA genes (n=%d)" , length(ahba_genes)),
#     color = "logFold" 
#   ) +
#   theme_bw(base_size = 14) +  # Clean layout
#   theme(
#     panel.grid.minor = element_blank(),
#     legend.position = c(0.80, 0.36),
#     legend.title = element_text(size = 12),
#     legend.key.height = unit(0.4, "cm"),
#     legend.key.width = unit(0.4, "cm"),
#     axis.title = element_blank(),
#     axis.text = element_text(size = 12),
#     plot.title = element_blank()  # Box around the plot
#   ) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#3f3f3f") +
#   # add p=0.05 line for x and y
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "#3f3f3f") +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "#3f3f3f") +
#   coord_fixed(1, xlim=c(0,1), ylim=c(0,1))
#   ggsave(sprintf('%s/Scatter_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.png', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test),
#             plot = p.psig, width = 3, height = 3, dpi = 300)
#  }
# }




# Create a list to store plots
plots <- list()
# Iterate through donors and percentage values
  for (k_perc in c(10, 20, 30)) {
    for (rdonor in c("0.6","0.4", "0.2", "-1")) {
  ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
    k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))
    # Determine if the plot is in the first column or last row
    show_x_axis <- (k_perc == 30)  # Show x-axis only for the last row (k_perc == 30)
    show_y_axis <- (rdonor == "0.6")  # Show y-axis only for the first column (rdonor == "0.2")
    show_title <- FALSE  # Show title only for the first column (rdonor == "0.2")
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
        title = if (show_title) sprintf("AHBA (n=%d)", length(ahba_genes)) else NULL
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
    plots[[paste(rdonor, k_perc)]] <- p.psig
  }
}

 final_p=wrap_plots(plots, nrow=3, ncol = 4) + 
 plot_annotation(title = sprintf("Background: AHBA vs. %s genes", anno_text),
  theme = theme(plot.title = element_text(size = 14, hjust = 0.5,face = "bold")))

ggsave(sprintf('%s/full_Scatter_plots_Anno_%s.png', res_path, anno2test),final_p, width = 12, height = 12, dpi = 300)


# Create a list to store plots
plots <- list()
# Iterate through donors and percentage values
  for (k_perc in c(10, 20)) {
    for (rdonor in c("0.6","0.4", "0.2")) {
  ahba_genes <- readRDS(sprintf('%s/AHBA_genes_rdonor%s.rds', data_path, rdonor))
    k <- round(k_perc * length(ahba_genes) / 100)
    res.sigfile <- sprintf('%s/Psig_sampled_perc%d_%d_rdnore%s_ahba%d_rep%d_Anno_%s.csv', res_path, k_perc, k, rdonor, length(ahba_genes), n_reps, anno2test)
    res.psig <- read.csv(res.sigfile) %>% mutate(ratio = (psig.protein+epsilon)/(psig.ahba+epsilon))

    # Determine if the plot is in the first column or last row
    show_x_axis <- (k_perc == 20)  # Show x-axis only for the last row (k_perc == 30)
    show_y_axis <- (rdonor == "0.6")  # Show y-axis only for the first column (rdonor == "0.2")
    show_title <- FALSE  # Show title only for the first column (rdonor == "0.2")
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
        title = if (show_title) sprintf("AHBA (n=%d)", length(ahba_genes)) else NULL
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
    plots[[paste(rdonor, k_perc)]] <- p.psig
  }
}

 final_p=wrap_plots(plots, nrow=2,ncol = 3) + 
 plot_annotation(title = sprintf("Background: AHBA vs. %s genes", anno_text),
  theme = theme(plot.title = element_text(size = 14, hjust = 0.5,face = "bold")))

 ggsave(sprintf('%s/selected_Scatter_plots_Anno_%s.png', res_path, anno2test),final_p, width = 9, height = 9, dpi = 300)
