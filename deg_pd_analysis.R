#!/usr/bin/env Rscript

# Load necessary libraries
library(edgeR)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "RNA-Seq Differential Expression Analysis")
parser$add_argument("--count-data-file", type = "character", default = "data_pTreg.csv", help = "Path to the count data file")
parser$add_argument("--output-prefix", type = "character", default = "vehicle_vs_evhd", help = "Prefix for output files")
parser$add_argument("--group-labels", type = "character", default = "Group1,Group2,Group3,Group4", help = "Comma-separated list of group labels")
parser$add_argument("--condition-labels", type = "character", default = "vehicle,evhd", help = "Comma-separated list of condition labels")
parser$add_argument("--n-replicates-per-group", type = "integer", default = 3, help = "Number of replicates per group")
parser$add_argument("--n-groups", type = "integer", default = 4, help = "Number of groups")
parser$add_argument("--ref-condition", type = "character", default = "vehicle", help = "Reference condition for comparison")
parser$add_argument("--target-condition", type = "character", default = "evhd", help = "Target condition for comparison")

# Parse arguments
args <- parser$parse_args()

# Function to perform differential expression analysis
run_deg_analysis <- function(count_data_file, output_prefix, group_labels, condition_labels, n_replicates_per_group, n_groups, ref_condition, target_condition) {
  
  # Expand group and condition labels based on the number of replicates
  group <- factor(rep(group_labels, each = n_replicates_per_group))
  condition <- factor(rep(condition_labels, n_groups))

  # Load the count data (ensure rownames are gene IDs)
  data <- read.csv(count_data_file, sep = ',', row.names = 'Geneid')
  
  # Filter to keep only the specified conditions
  selected <- condition %in% c(ref_condition, target_condition)
  data <- data[, selected]
  group <- group[selected]
  condition <- condition[selected]

  # Relevel the condition to set the reference
  condition <- relevel(condition, ref = ref_condition)
  
  # Create the design matrix for differential expression analysis
  design <- model.matrix(~ group + condition)
  
  # Create a DGEList object from the counts
  dge <- DGEList(counts = data)
  
  # Normalize counts using TMM normalization
  dge <- calcNormFactors(dge)

  # Estimate dispersion
  dge <- estimateDisp(dge, design)

  # Fit the generalized linear model
  fit <- glmFit(dge, design)

  # Perform the likelihood ratio test for the comparison
  lrt <- glmLRT(fit, coef = paste0("condition", target_condition))

  # Extract DEGs with all statistics
  de_genes <- topTags(lrt, n = Inf)$table
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  de_genes$pvals_adj <- p.adjust(de_genes$PValue, method = "BH")
  
  # Save DEGs to a CSV file
  write.csv(de_genes, file = paste0(output_prefix, "_DEGs.csv"), quote = TRUE, row.names = TRUE)
  
  return(de_genes)
}

# Run the analysis with command-line arguments
deg_results <- run_deg_analysis(
  count_data_file = args$count_data_file,
  output_prefix = args$output_prefix,
  group_labels = strsplit(args$group_labels, ",")[[1]],  # Convert comma-separated string to vector
  condition_labels = strsplit(args$condition_labels, ",")[[1]],  # Convert comma-separated string to vector
  n_replicates_per_group = args$n_replicates_per_group,
  n_groups = args$n_groups,
  ref_condition = args$ref_condition,
  target_condition = args$target_condition
)
