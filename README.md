# Synteny Plot Script

This R script creates a circular synteny plot from your final synteny CSV files generated from the [synteny_finder python script](https://github.com/ian-bda/synteny_finder), automatically detecting and highlighting CD300 genes.

<img width="787" height="792" alt="Screenshot 2025-08-12 at 10 39 47 AM" src="https://github.com/user-attachments/assets/1c0aefca-99d7-408d-881b-119dd00aa7ae" />



## Features

- **Dynamic Data Loading**: Automatically reads all CSV files from the `final_synteny_csvs` directory
- **CD300 Gene Detection**: Automatically identifies and highlights CD300 genes (case-insensitive)
- **Fallback Highlighting**: If no CD300 genes are found, highlights the top-scoring genes
- **Flexible Configuration**: Easy-to-modify parameters at the top of the script
- **PDF Output**: Saves the plot as a high-quality PDF file
- **Comprehensive Summary**: Provides detailed statistics about the analysis

## Requirements

- R with the following packages:
  - `circlize`
  - `dplyr`
  - `readr`

## Installation

```r
install.packages(c("circlize", "dplyr", "readr"))
```

## Usage

1. Place your synteny CSV files in the `final_synteny_csvs` directory
2. Run the script:
   ```bash
   Rscript Macrosynteny_plot.R
   ```
3. The plot will be saved as `synteny_plot.pdf` in the current directory

## Configuration

Edit the parameters at the top of the script to customize:

```r
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
```

## Input File Format

The script expects tab-separated CSV files with the following columns:
- `gene1`: Gene name from species 1
- `gene2`: Gene name from species 2
- `score`: Synteny score
- `chrom1`: Chromosome from species 1
- `start1`: Start position from species 1
- `end1`: End position from species 1
- `chrom2`: Chromosome from species 2
- `start2`: Start position from species 2
- `end2`: End position from species 2

## Output

- **PDF Plot**: Circular synteny visualization with highlighted genes
- **Console Summary**: Statistics including:
  - Total synteny links
  - Number of highlighted links
  - Number of background links
  - Species pairs analyzed
  - Chromosomes included
  - Genes highlighted
  - CD300 genes found (if any)

## Example Output

```
=== SYNTHENY PLOT SUMMARY ===
Total synteny links: 30168 
Highlighted links: 214 
Background links: 29954 
Species pairs analyzed: 45 
Chromosomes included: 8 
Genes highlighted: 32 
CD300 genes found: 32 
```

## Troubleshooting

- **Empty files**: The script automatically skips empty CSV files
- **Missing columns**: Files with incorrect column structure are skipped
- **No data**: The script will stop with an error message if no valid data is found
- **Plot warnings**: Some "out of plotting region" warnings are normal and don't affect the final plot

## Customization

To modify the gene highlighting pattern, edit the `identify_cd300_genes` function:

```r
identify_cd300_genes <- function(genes) {
  # Change this pattern to highlight different genes
  cd300_pattern <- "cd300[a-z]*"
  cd300_genes <- genes[grepl(cd300_pattern, genes, ignore.case = TRUE)]
  return(cd300_genes)
}
```
