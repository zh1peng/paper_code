library(viridis)
library(scales)
library(ggplot2)
library(gghalves)
library(dplyr)
# install.packages('stringdist')
library(stringdist)
library(dplyr)


data_path = 'E:/xhmhc/ORA misuse/data'
res_path = 'E:/xhmhc/ORA misuse/article_res'

if (!dir.exists(res_path)) {
  dir.create(res_path)
}

df = readxl::read_excel(file.path(data_path, 'extracted_information.xlsx'))

year_counts <- df %>%
  group_by(Year) %>%
  summarise(Count = n()) %>% 
  ungroup()

# color2use = viridis(4)[4]
# Create the plot
p <- ggplot(year_counts, aes(x = Year, y = Count)) +
  geom_line(color ='#fde725' , size = 2) + # Use a harmonized color for the line (from viridis palette)
  geom_point(color = '#3c9390', size = 3) + # Use a contrasting harmonized color for the points
  theme_minimal(base_size = 15) +
  labs(title = "Trend of Research Articles Over Years",
       x = "Year",
       y = "Number of Articles") +
  theme(
    plot.title = element_blank(), # Optionally blank out title
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "grey80", linetype = "dashed"), # Show only major gridlines
    panel.grid.minor = element_blank(), # Hide minor gridlines
    axis.text = element_text(size = 12,color = "black"),
  ) +
  scale_x_continuous(breaks = seq(min(year_counts$Year), max(year_counts$Year), 1)) +
  scale_y_continuous(labels = comma)

# Save the plot with specified dimensions and resolution
ggsave(file.path(res_path, 'year_counts_harmonized.png'), plot = p, width = 4, height = 3, dpi = 300)



# plot journal counts
journal_counts <- df %>%
  group_by(Journal) %>%
  summarise(Count = n()) %>%
  ungroup()

# Remove "(New York, N.Y.:1991)" from "Cerebral Cortex" for better display
journal_counts$Journal <- gsub("\\s*\\(New York, N\\.Y\\. : 1991\\)", "", journal_counts$Journal)

# Get the top 5 journals
top_5_journals <- journal_counts[order(-journal_counts$Count), ][1:5, ]

# Calculate the percentage for each of the top 5 journals
top_5_journals$Percentage <- (top_5_journals$Count / sum(journal_counts$Count)) * 100

# Create the horizontal bar plot
p=ggplot(top_5_journals, aes(x = Count, y = reorder(Journal, Count))) +
  geom_bar(stat = "identity", fill = "#fde725") + 
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), hjust = 1, size = 4) +
  scale_y_discrete(position = "right") +
  labs(title = "Top 5 Journals by Percentage of Articles",
       x = "Number of Articles",
       y = "Journal") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_text(size =14),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    axis.text = element_text(size = 12,color = "black"),
  ) 

ggsave(file.path(res_path, 'top_5_journals.png'), plot = p, width = 3.5, height = 2.75, dpi = 300)


# make plot for tools
# Extract the tools column and handle multiple tools
tools_split <- stringr::str_split(df$Tool, pattern = "/", simplify = FALSE)

# Unlist the tools into a single vector and create a data frame
tools_vector <- unlist(tools_split)
tools_counts <- as.data.frame(table(tools_vector))
colnames(tools_counts) <- c("Tool", "Count")

# Get the top 5 tools
top_tools <- tools_counts[order(-tools_counts$Count), ][1:8, ]

# Calculate percentages for the top 5 tools
top_tools$Percentage <- (top_tools$Count / sum(tools_counts$Count)) * 100

# Create the horizontal bar plot
p <- ggplot(top_tools, aes(x = Count, y = reorder(Tool, Count))) +
  geom_bar(stat = "identity", fill = "#fde725") + 
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), hjust = 1, size = 4) +
  labs(title = "Top 5 Tools by Usage Percentage",
       x = "Count",
       y = "Tool") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12,color = "black"),
  )

# Save the plot
ggsave(file.path(res_path, 'top_tools.png'), plot = p, width = 4, height = 3, dpi = 300)



## plot bg genes used for analysis
# Convert the 'N genes' column to numeric where possible
df$N_genes_numeric <- as.numeric(df$`N genes`) # Converts numeric values, others become NA

# Count "Not explicitly stated" occurrences
not_stated_count <- sum(is.na(df$N_genes_numeric))
# Load necessary libraries




# Filter the numeric data (exclude NA values)
numeric_data <- df %>%
  mutate(N_genes_numeric = as.numeric(`N genes`)) %>%
  filter(!is.na(N_genes_numeric))

median(numeric_data$N_genes_numeric)
mean(numeric_data$N_genes_numeric)

# Create the raincloud plot
p <- ggplot(numeric_data, aes(x = "", y = N_genes_numeric)) +
  geom_half_violin(aes(fill = "Density"), side = "r", alpha = 0.7, color = NA) + # Density plot (half-violin)
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black", alpha = 0.5) + # Box plot
  geom_jitter(aes(color = "Points"), width = 0.1, alpha = 0.7, size = 1.5) + # Scatter plot (jittered points)
  scale_fill_manual(values = c("Density" = "#21908C")) +
  scale_color_manual(values = c("Points" = "#440154")) +
  labs(title = "Raincloud Plot for N Genes",
       x = NULL,
       y = "N Genes") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Save the plot
ggsave(file.path(res_path, 'raincloud_plot_n_genes.png'), plot = p, width = 8, height = 6, dpi = 300)






numeric_data <- df %>%
  mutate(N_genes_numeric = as.numeric(`N genes`)) %>% # Convert N_genes to numeric
  filter(!is.na(N_genes_numeric)) %>%# Filter out rows with NA values
  mutate(group = case_when(
    report %in% c("state_but_different", "state_but_less", "state_but_more") & !is.na(N_genes_numeric) ~ "Reported but different",
    report == "stated_same" ~ "Reported same",
    report == "not_stated" ~ "Not reported",
    TRUE ~ NA_character_ # Exclude rows that don't match the criteria
  )) %>%
  filter(!is.na(group)) # Remove rows without valid group classification



col2use  = viridis(4)[c(4,2,3)] 
# Create the raincloud plot
p <- ggplot(numeric_data, aes(x = group, y = N_genes_numeric, fill = group)) +
  # geom_half_violin(side = "l", alpha = 0.7, color = NA) + # Density (half violin)
  geom_point(aes(color = group), position = position_jitterdodge(), size = 1, alpha = 0.8) + # Jittered points
  geom_boxplot(width = 0.15, outlier.shape = NA,  color = "black", alpha = 0.6) + # Boxplot
  scale_fill_manual(values=col2use) +
  scale_color_manual(values=col2use) +
  labs(title = "Raincloud Plot of N Genes by Group",
       x = "Group",
       y = "AHBA Genes") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 14),
    axis.text.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12,color = "black")
  ) 
ggsave(file.path(res_path, 'raincloud_plot_n_genes_by_group.png'), plot = p, width = 2.2, height =1.5, dpi = 300)


# create raincloud plot for sig %

numeric_data <- df %>%
  mutate(Sig_genes_numeric = as.numeric(`Significant genes percentage`)) %>% # Convert N_genes to numeric
  filter(!is.na(Sig_genes_numeric)) %>%# Filter out rows with NA values
  mutate(group = case_when(
    report %in% c("state_but_different", "state_but_less", "state_but_more") & !is.na(Sig_genes_numeric) ~ "Reported but different",
    report == "stated_same" ~ "Reported same",
    report == "not_stated" ~ "Not reported",
    TRUE ~ NA_character_ # Exclude rows that don't match the criteria
  )) %>%
  filter(!is.na(group)) # Remove rows without valid group classification

col2use  = viridis(4)[c(4,2,3)] 
# Create the raincloud plot
p <- ggplot(numeric_data, aes(x = group, y = Sig_genes_numeric, fill = group)) +
  # geom_half_violin(side = "l", alpha = 0.7, color = NA) + # Density (half violin)
  geom_point(aes(color = group), position = position_jitterdodge(), size = 1, alpha = 0.8) + # Jittered points
  geom_boxplot(width = 0.15, outlier.shape = NA,  color = "black", alpha = 0.6) + # Boxplot
  scale_fill_manual(values=col2use) +
  scale_color_manual(values=col2use) +
  labs(title = "Raincloud Plot of N Genes by Group",
       x = "Group",
       y = "Sig Genes (%)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 14),
    axis.text.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12,color = "black")
  ) 
ggsave(file.path(res_path, 'raincloud_plot_sig_genes_by_group.png'), plot = p, width = 2, height =1.5, dpi = 300)


# Create the pie chart
# Summarize the count of each group
numeric_data <- df %>%
  mutate(group = case_when(
    report %in% c("state_but_different", "state_but_less", "state_but_more")  ~ "Reported but different",
    report == "state_but_no_ahba" ~ "Reported but no AHBA",
    report == "stated_same" ~ "Reported same",
    report == "not_stated" ~ "Not reported",
    TRUE ~ NA_character_ # Exclude rows that don't match the criteria
  )) 

group_counts <- numeric_data %>%
  group_by(group) %>%
  summarize(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) # Calculate percentages

col2use  = viridis(4)[c(4,2,1,3)] 
pie_chart <- ggplot(group_counts, aes(x = "", y = Percentage, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y",  start = pi*0.4) + # Convert to pie chart
  scale_fill_manual(values=col2use) + # Use viridis colors
  labs(title = "Distribution of Correct Categories",
       y = NULL,
       x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  )
  # ) +
  #  geom_text_repel(
  #   aes(label = paste0(group, ": ", sprintf("%.1f%%", Percentage))),
  #   #position = position_stack(vjust = 0.5),
  #   nudge_x = 0.8, # Nudge labels outside the pie
  #   nudge_y = -1,
  #   direction = "y",
  #   box.padding = 1,
  #   size = 4,
  #   segment.color = "black" # Line connecting label and pie
  # )
# Save the pie chart
ggsave(file.path(res_path, 'pie_chart_correct_categories.png'), plot = pie_chart, width = 6, height = 6, dpi = 300)



#
terms_split <- stringr::str_split(df$`key enrichment`, pattern = ",", simplify = FALSE)

# Unlist the tools into a single vector and create a data frame
terms_vector <- unlist(terms_split)

# strip leading and trailing white spaces
terms_vector <- trimws(terms_vector)
terms_vector <- tolower(terms_vector)



data <- as.data.frame(table(terms_vector))
colnames(data) <- c("term", "count")

# write to csv to do manual curation
# write.csv(data, file.path(res_path, 'terms_count.csv'), row.names = FALSE)

new_terms =readxl::read_excel(file.path(data_path, 'terms.xlsx'))
data$new_terms <- new_terms$Terms
data <- data[-1,] 

lookup <- setNames(data$new_terms, data$term)
terms_vector <- replace(terms_vector, terms_vector %in% names(lookup), lookup[terms_vector])

data <- as.data.frame(table(terms_vector))
colnames(data) <- c("term", "count")
# write.csv(data, file.path(res_path, 'terms_count2.csv'), row.names = FALSE)

# # Compute pairwise string distances
# dist_matrix <- stringdistmatrix(data$term, data$term, method = "jw")  # Jaro-Winkler distance
# # Perform hierarchical clustering
# hc <- hclust(as.dist(dist_matrix), method = "average")
# # Cut tree into groups based on a similarity threshold
# groups <- cutree(hc, h = 0.3)  # Adjust threshold (e.g., 0.2 for similar terms)
# # Add group labels to the original data
# data$group <- groups
# # Summarize counts by group
# grouped_data <- data %>%
#   group_by(group) %>%
#   summarise(
#     combined_terms = paste(unique(term), collapse = "; "),
#     total_count = sum(count)
#   ) %>%
#   arrange(desc(total_count))
# # View grouped data
# print(grouped_data)


library(wordcloud)

wordcloud(
  words = data$term,            # Words to display
  freq = data$count,            # Frequencies of words
  min.freq = 5,                 # Minimum frequency to include a word
  max.words = 1000,              # Maximum number of words to include
  random.order = FALSE,         # Display most frequent words first
  rot.per = 0.2,               # Proportion of words with vertical rotation
  colors = viridis(100) # Color palette
)

data2save <- data %>% filter(count > 1) %>%
  arrange(desc(count)) %>%
  rename(weight = count, word = term) %>%
  mutate(color = substr(                  # Extract standard hex code (first 7 characters)
      viridis(20)[cut(weight, breaks = 20)],  # Generate color based on weight
      1, 7                          # Take the first 7 characters of the hex
    )
  ) %>%
  select(weight, word, color)
write.csv(data2save, file.path(res_path, 'wordcloud_data.csv'), row.names = FALSE)




