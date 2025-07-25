library(circlize)
library(dplyr)
library(readr)

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Set the directory containing your CSV files
DATA_DIR <- '/Users/ianbda/Documents/NC_State_Genomics/Yoder_Lab/CD300_Project/Euarchontoglires/synteny/final_synteny_csvs'

# Output settings
SAVE_PLOT <- TRUE
OUTPUT_FILE <- 'synteny_plot.pdf'
PLOT_WIDTH <- 12
PLOT_HEIGHT <- 12

# Plot appearance settings
BACKGROUND_COLOR <- "#80808020"  # Semi-transparent gray
HIGHLIGHT_COLOR <- "#FF0000"     # Red for highlighted genes
BACKGROUND_LINE_WIDTH <- 0.5
HIGHLIGHT_LINE_WIDTH <- 3
GAP_DEGREE <- 5

# =============================================================================
# MAIN SCRIPT
# =============================================================================

# Set working directory
setwd(DATA_DIR)

# List all CSV files
file_list <- list.files(pattern = "*.csv")

# Function to extract species names from filename
extract_species_from_filename <- function(filename) {
  # Remove "synteny_" prefix and "_parsed.csv" suffix
  clean_name <- gsub("^synteny_", "", filename)
  clean_name <- gsub("_parsed\\.csv$", "", clean_name)
  
  # Split by "." to get the two species
  species_parts <- strsplit(clean_name, "\\.")[[1]]
  
  if (length(species_parts) == 2) {
    # Clean up species names by removing "_parsed" suffix if present
    species1 <- gsub("_parsed$", "", species_parts[1])
    species2 <- gsub("_parsed$", "", species_parts[2])
    
    return(list(species1 = species1, species2 = species2))
  } else {
    return(NULL)
  }
}

# Read and combine all synteny data
synteny_data <- data.frame()

for (file in file_list) {
  # Skip empty files
  file_content <- readLines(file)
  if (length(file_content) <= 1) {
    next
  }
  
  # Read CSV with proper column types (tab-separated)
  temp_data <- read_tsv(file, col_types = cols(
    gene1 = col_character(),
    gene2 = col_character(),
    score = col_double(),
    chrom1 = col_character(),
    start1 = col_double(),
    end1 = col_double(),
    chrom2 = col_character(),
    start2 = col_double(),
    end2 = col_double()
  ))
  
  # Check if the data frame has the expected structure
  if (!all(c("gene1", "gene2", "chrom1", "chrom2", "start1", "start2") %in% names(temp_data))) {
    next
  }
  
  # Add species information
  species_info <- extract_species_from_filename(file)
  if (!is.null(species_info)) {
    temp_data$species1 <- species_info$species1
    temp_data$species2 <- species_info$species2
    
    # Normalize chromosome names by removing ".0" suffixes
    temp_data$chrom1_clean <- gsub("\\.0$", "", temp_data$chrom1)
    temp_data$chrom2_clean <- gsub("\\.0$", "", temp_data$chrom2)
    
    # Create unique chromosome identifiers that include species information
    temp_data$chrom1_unique <- paste0(temp_data$species1, "_chr", temp_data$chrom1_clean)
    temp_data$chrom2_unique <- paste0(temp_data$species2, "_chr", temp_data$chrom2_clean)
  }
  
  # Combine with main dataset
  if (nrow(synteny_data) == 0) {
    synteny_data <- temp_data
  } else {
    synteny_data <- rbind(synteny_data, temp_data)
  }
}

# Check if we have any data
if (nrow(synteny_data) == 0) {
  stop("No synteny data found in any CSV files!")
}

# Remove any rows with missing data
synteny_data <- synteny_data %>%
  filter(!is.na(chrom1) & !is.na(chrom2) & !is.na(start1) & !is.na(start2))

# Check if we still have data after filtering
if (nrow(synteny_data) == 0) {
  stop("No valid synteny data found after filtering!")
}

# Get all unique chromosomes dynamically (now using unique chromosome identifiers)
all_chromosomes <- unique(c(synteny_data$chrom1_unique, synteny_data$chrom2_unique))

# Get all unique genes for potential highlighting
all_genes <- unique(c(synteny_data$gene1, synteny_data$gene2))

# Function to identify CD300 genes (case-insensitive)
identify_cd300_genes <- function(genes) {
  cd300_pattern <- "cd300[a-z]*"
  cd300_genes <- genes[grepl(cd300_pattern, genes, ignore.case = TRUE)]
  return(cd300_genes)
}

# Find CD300 genes automatically
cd300_genes <- identify_cd300_genes(all_genes)

# If no CD300 genes found, use genes with highest scores as highlights
if (length(cd300_genes) == 0) {
  # Get top scoring genes for highlighting
  top_genes <- synteny_data %>%
    group_by(gene1, gene2) %>%
    summarise(max_score = max(score, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(max_score)) %>%
    head(20)
  
  highlight_genes <- unique(c(top_genes$gene1, top_genes$gene2))
} else {
  highlight_genes <- cd300_genes
}

# Scale positions to reasonable values (convert to Mb)
synteny_data$start1_mb <- synteny_data$start1 / 1e6
synteny_data$start2_mb <- synteny_data$start2 / 1e6

# Calculate chromosome ranges dynamically (using unique chromosome identifiers)
chrom_ranges <- data.frame(
  chrom = all_chromosomes,
  start = 0,
  end = sapply(all_chromosomes, function(chr) {
    max(c(synteny_data$start1_mb[synteny_data$chrom1_unique == chr],
          synteny_data$start2_mb[synteny_data$chrom2_unique == chr]), na.rm = TRUE)
  })
)

# Ensure minimum chromosome size for visualization
min_chrom_size <- 1.0  # 1 Mb minimum
chrom_ranges$end <- pmax(chrom_ranges$end, min_chrom_size)

# Normalize chromosome lengths for circular plot
chrom_ranges$length <- chrom_ranges$end - chrom_ranges$start
total_length <- sum(chrom_ranges$length)
chrom_ranges$prop <- chrom_ranges$length / total_length
chrom_ranges$end_norm <- cumsum(chrom_ranges$prop) * 100
chrom_ranges$start_norm <- c(0, chrom_ranges$end_norm[-nrow(chrom_ranges)])

# Create the xlim matrix for circos.initialize
xlim_matrix <- cbind(chrom_ranges$start_norm, chrom_ranges$end_norm)
rownames(xlim_matrix) <- chrom_ranges$chrom

# Start PDF device if saving is enabled
if (SAVE_PLOT) {
  pdf(OUTPUT_FILE, width = PLOT_WIDTH, height = PLOT_HEIGHT)
}

# Initialize the circular layout
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.degree = GAP_DEGREE)

circos.initialize(factors = chrom_ranges$chrom, xlim = xlim_matrix)

# Add base ideogram track
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05)

# Add chromosome labels outside the sections (with species information)
circos.track(track.index = 1, panel.fun = function(x, y) {
  # Extract species and chromosome info for display
  chrom_label <- CELL_META$sector.index
  # Remove the "_chr" part and format nicely
  display_label <- gsub("_chr", " chr", chrom_label)
  circos.text(CELL_META$xcenter,
              CELL_META$ylim[2] + mm_y(10),
              display_label,
              cex = 0.7,
              niceFacing = TRUE)
}, bg.border = NA)

# Separate highlighted and non-highlighted links
highlighted_links <- synteny_data %>%
  filter(gene1 %in% highlight_genes | gene2 %in% highlight_genes)

non_highlighted_links <- synteny_data %>%
  filter(!(gene1 %in% highlight_genes | gene2 %in% highlight_genes))

# Function to normalize position within chromosome (updated for unique chromosomes)
normalize_position <- function(pos, chrom, xlim_matrix, data, chrom_col, pos_col) {
  chr_data <- data[data[[chrom_col]] == chrom, ]
  if (nrow(chr_data) == 0) return(0)
  
  max_pos <- max(chr_data[[pos_col]], na.rm = TRUE)
  if (max_pos == 0) return(0)
  
  start_norm <- xlim_matrix[chrom, 1]
  end_norm <- xlim_matrix[chrom, 2]
  
  normalized <- (pos / max_pos) * (end_norm - start_norm) + start_norm
  return(normalized)
}

# Plot non-highlighted links first (background)
if (nrow(non_highlighted_links) > 0) {
  for (i in seq_len(nrow(non_highlighted_links))) {
    chr1 <- as.character(non_highlighted_links$chrom1_unique[i])
    chr2 <- as.character(non_highlighted_links$chrom2_unique[i])
    
    start1_norm <- normalize_position(
      non_highlighted_links$start1_mb[i], 
      chr1, 
      xlim_matrix, 
      synteny_data, 
      "chrom1_unique", 
      "start1_mb"
    )
    
    start2_norm <- normalize_position(
      non_highlighted_links$start2_mb[i], 
      chr2, 
      xlim_matrix, 
      synteny_data, 
      "chrom2_unique", 
      "start2_mb"
    )
    
    circos.link(chr1, start1_norm, chr2, start2_norm,
                col = BACKGROUND_COLOR, lwd = BACKGROUND_LINE_WIDTH)
  }
}

# Plot highlighted links on top
if (nrow(highlighted_links) > 0) {
  for (i in seq_len(nrow(highlighted_links))) {
    chr1 <- as.character(highlighted_links$chrom1_unique[i])
    chr2 <- as.character(highlighted_links$chrom2_unique[i])
    
    start1_norm <- normalize_position(
      highlighted_links$start1_mb[i], 
      chr1, 
      xlim_matrix, 
      synteny_data, 
      "chrom1_unique", 
      "start1_mb"
    )
    
    start2_norm <- normalize_position(
      highlighted_links$start2_mb[i], 
      chr2, 
      xlim_matrix, 
      synteny_data, 
      "chrom2_unique", 
      "start2_mb"
    )
    
    circos.link(chr1, start1_norm, chr2, start2_norm,
                col = HIGHLIGHT_COLOR, lwd = HIGHLIGHT_LINE_WIDTH)
  }
}

# Add ticks and labels
for (chr in all_chromosomes) {
  start_pos <- xlim_matrix[chr, 1]
  end_pos <- xlim_matrix[chr, 2]
  
  tick_positions <- seq(start_pos, end_pos, length.out = 5)
  labels <- c("0", "25", "50", "75", "100")
  
  circos.axis(
    sector.index = chr,
    major.at = tick_positions,
    labels = labels,
    labels.cex = 0.6,
    labels.facing = "outside"
  )
}

# Close PDF device if saving is enabled
if (SAVE_PLOT) {
  dev.off()
  cat("Plot saved to:", OUTPUT_FILE, "\n")
}

# Print summary statistics
cat("\n=== SYNTHENY PLOT SUMMARY ===\n")
cat("Total synteny links:", nrow(synteny_data), "\n")
cat("Highlighted links:", nrow(highlighted_links), "\n")
cat("Background links:", nrow(non_highlighted_links), "\n")
cat("Species pairs analyzed:", length(file_list), "\n")
cat("Unique chromosomes included:", length(all_chromosomes), "\n")
cat("Genes highlighted:", length(highlight_genes), "\n")
if (length(cd300_genes) > 0) {
  cat("CD300 genes found:", length(cd300_genes), "\n")
}

# Print unique chromosomes for reference
cat("\nUnique chromosomes in the plot:\n")
for (chr in sort(all_chromosomes)) {
  display_name <- gsub("_chr", " chr", chr)
  cat("  ", display_name, "\n")
}

