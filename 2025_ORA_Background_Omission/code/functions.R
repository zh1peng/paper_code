
read_AHBA_data <- function(folder_path='E:/xhmhc/ORA misuse/data/AHBA_ENIGMA', 
                           rdonor = c("-1", "0.2", "0.4", "0.6", "0.8"), 
                           atlas = c("desikan", "glasser_360", 
                                     "schaefer_100", "schaefer_200", 
                                     "schaefer_300", "schaefer_400")) {
  # Validate the atlas and rdonor input
  atlas <- match.arg(atlas)
  rdonor <- match.arg(rdonor)
  
  # Construct the search pattern for atlas
  atlas_pattern <- switch(
    atlas,
    "desikan" = "",  # Matches anything for desikan
    atlas # Exact match for other atlases
  )
  # Construct the search pattern for rdonor
  rdonor_pattern <- paste0("r", rdonor)
 
  # Build the full pattern to search for a matching file
  csv_file <- paste0("allgenes_stable_", rdonor_pattern, atlas_pattern, ".csv")
  data <- read.csv(file.path(folder_path,csv_file))
  return(data) # Return the data frame
}

fast_ORA <- function(sig_genes, gs_list, bg_genes) {
  # Validate inputs
  if (!is.character(sig_genes) || !is.list(gs_list) || !is.character(bg_genes)) {
    stop("Invalid input types: 'sig_genes' and 'bg_genes' must be character vectors, and 'gs_list' must be a named list of character vectors.")
  }
  
  if (length(sig_genes) == 0 || length(gs_list) == 0 || length(bg_genes) == 0) {
    stop("Inputs cannot be empty. Please provide non-empty 'sig_genes', 'gs_list', and 'bg_genes'.")
  }
  
  if (is.null(names(gs_list))) {
    stop("'gs_list' must be a named list with descriptive names for each gene set.")
  }
  
  # Intersect sig_genes and gene sets with the background
  sig_genes <- intersect(sig_genes, bg_genes)
  gs_list <- intersectToList(gs_list, bg_genes)
  
  # Get background size and query size
  n_bg <- length(bg_genes)
  n_sig <- length(sig_genes)
  
  # Calculate overlap, gene set sizes, and background sizes
  x <- sapply(intersectToList(gs_list, sig_genes), length) # Overlap sizes
  m <- sapply(gs_list, length)                            # Gene set sizes
  n <- n_bg - m                                           # Background sizes
  k <- n_sig                                              # Query size
  
  # Perform one-sided hypergeometric test
  p_values <- phyper(x - 1, m, n, k, lower.tail = FALSE)
  names(p_values) <- names(gs_list) # Retain gene set names
  # Return p-values
  return(p_values)
}


# Parallel simulation function
simulate_repetition <- function(rep_id, k, selected.gs, ahba_genes,protein_genes) {
  # Simulate significant genes
  sig_genes <- sample(ahba_genes, k)
  # Compute p-values for the two analyses
  pvals.protein <- fast_ORA(sig_genes, selected.gs, protein_genes)
  pvals.ahba <- fast_ORA(sig_genes, selected.gs, ahba_genes)
  # Return results as a data frame
  data.frame(
    Pathway = names(selected.gs),
    pvals.protein = pvals.protein,
    pvals.ahba = pvals.ahba,
    Rep = rep_id  # Track the repetition ID
  )
}

