library(rmarkdown)

# Set the working directory to the directory containing the R Markdown files
setwd("F:/Google Drive/post-doc/vitural_histology_revisit/revision_code/analysis")

# Get a list of all R Markdown files in the directory
rmd_files <- list.files(pattern = "\\.Rmd$")

# Render each R Markdown file
for (file in rmd_files) {
  # Render the R Markdown file
  render(input = file)
}