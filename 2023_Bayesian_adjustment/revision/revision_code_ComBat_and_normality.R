library(dplyr)
library(tidyr)
library(ggplot2)
data_path='F:/Google Drive/post-doc/Bayesian_Project/new_model/revision'

# Initialize an empty list to store p-values
site_p <- list()
site_t <- list()
site_d <- list()
site_diff_mean <- list()
site_diff_sd <- list()
site_minial_sd <- list()
# Assuming data_path is predefined somewhere in your code
# data_path <- "your/data/path"

for (site_i in 1:32) {
  # Try to read the files; if there's an error, just print the message
  tryCatch({
    sampled_es_df_combat <- read.csv(sprintf('%s/site_data/site_%s_sampled_es.csv', data_path, site_i))
    sampled_es_df_nocombat <- read.csv(sprintf('%s/site_data_no_combat/site_%s_sampled_es.csv', data_path, site_i))
    
    # Perform the operations only if both files are read successfully
    site_var_combat <- apply(sampled_es_df_combat[, grep('L_|R_', names(sampled_es_df_combat))], 2, sd)
    site_var_no_combat <- apply(sampled_es_df_nocombat[, grep('L_|R_', names(sampled_es_df_nocombat))], 2, sd)
    
    # Perform paired t-test and store the p-value
    tmp=t.test(site_var_combat, site_var_no_combat, paired = TRUE)
    site_p[[paste0('site_',site_i)]] <- tmp$p.value
    site_t[[paste0('site_',site_i)]] <- tmp$statistic

    # Calculate means of the differences
    mean_diff = mean(site_var_combat - site_var_no_combat)
    sd_diff = sd(site_var_combat - site_var_no_combat)
    cohen_d = mean_diff / sd_diff
    site_diff_mean[[paste0('site_',site_i)]] <- mean_diff
    site_diff_sd[[paste0('site_',site_i)]] <- sd_diff
    site_d[[paste0('site_',site_i)]] <- cohen_d
    site_minial_sd[[paste0('site_',site_i)]] <- mean(site_var_combat)
      }, error = function(e) {
        cat("An error occurred for site", site_i, ":", e$message, "\n")
        # No need to use 'next' here; just end the error function
      })
    }
site_n=as.numeric(readRDS('F:/Google Drive/post-doc/Bayesian_Project/new_model/site_N.rds')) #
site_p_df <- data.frame(site = names(site_p),
                        t_value = unlist(site_t),
                        p_value = unlist(site_p),
                        d_value = unlist(site_d),
                        site_diff_mean = unlist(site_diff_mean),
                        site_diff_sd = unlist(site_diff_sd),
                        site_minial_sd = unlist(site_minial_sd),
                        site_n=site_n)

# report site_minial_sd of the site_var_combat
mean(site_p_df$site_minial_sd)
sd(site_p_df$site_minial_sd)

# here we should look at the mean of the differences
# because that is directly related to the model

n <- 150  # sample size
site_p_df=site_p_df %>% 
mutate(site_diff_se = site_diff_sd / sqrt(n)) %>%
arrange(site_n) %>%
mutate(idx=1:n(),
       site_n2show=paste0('Site',idx,'(',site_n,')')) 

p=ggplot(site_p_df, aes(x=site_n2show, y=site_diff_mean)) + 
  geom_point(size=3) + 
  geom_errorbar(aes(ymin=site_diff_mean - site_diff_se, ymax=site_diff_mean + site_diff_se), width=.5) +
  theme_minimal()+
  xlab('Site (sample size)')+
  ylab('Mean of the differences')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-tick labels

ggsave('F:/Google Drive/post-doc/Bayesian_Project/new_model/revision/figures/mean_diff.png',p, 
      bg='white', width=10, height=5, units='in', dpi=300)



# test normality of the effec sizes
site_norm.p <- list()
# test normality of the data
for (site_i in 1:32) {
tryCatch({
    sampled_es_df_combat <- read.csv(sprintf('%s/site_data/site_%s_sampled_es.csv', data_path, site_i))
    site_norm.p[[paste0('site_',site_i)]]=apply(sampled_es_df_combat[, grep('L_|R_', names(sampled_es_df_combat))],2, function(data_vector) ks.test(data_vector, "pnorm", mean=mean(data_vector), sd=sd(data_vector))$p.value)
 

}, error = function(e) {
        cat("An error occurred for site", site_i, ":", e$message, "\n")
        # No need to use 'next' here; just end the error function
      })
}
 
df_site_norm_p_all <- do.call(rbind, lapply(site_norm.p, function(x) as.data.frame(t(x)))) %>%
                      mutate(site_n=site_n) %>%
                      arrange(site_n) %>%
                      mutate(idx=1:n(),
                             site_n2show=paste0('Site',idx,'(',site_n,')'))%>%
                      select(-idx,-site_n)%>%
                      pivot_longer(cols = -c(site_n2show), names_to = "ROI", values_to = "p_value")
# make violin plot for p_values for each site_n2show
p=ggplot(df_site_norm_p_all, aes(x=site_n2show, y=p_value)) + 
  geom_violin() + 
  geom_boxplot(width=0.1)+
  theme_minimal()+
  xlab('Site (sample size)')+
  ylab('p-value')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-tick labels




library(ggplot2)

# Define a vector of 21 distinct colors
colors <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", 
            "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
            "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", 
            "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
            "#ffffff", "#000000")

# Enhanced ggplot code with a custom color palette
ggplot(df_site_norm_p_all, aes(x=site_n2show, y=p_value, fill=site_n2show)) +
  geom_violin() + 
  geom_boxplot(width=0.1, outlier.shape = NA) +  # Hide outliers in boxplot to avoid duplication with jitter
  geom_jitter(width=0.2, size=1.5, shape=21, color="black", fill="white") +  # Add jittered points
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") +  # Horizontal line at p=0.05
  scale_fill_manual(values=colors) +  # Custom colors
  theme_minimal() +
  xlab('Site (sample size)') +
  ylab('p-value') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Hide legend if not needed
  ) 
ggsave('F:/Google Drive/post-doc/Bayesian_Project/new_model/revision/figures/p_value.png',p, 
      bg='white', width=10, height=5, units='in', dpi=300)




# Now each entry is a row in 'df_site_norm_p'


