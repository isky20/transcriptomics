# Load required libraries
library(edgeR)  # For differential expression analysis
library(dplyr)  # For data manipulation
library(argparse)  # For parsing command-line arguments

# Function to get Differentially Expressed Genes (DEGs)
get_DEG <- function(raw.reads.csv, colname.case, number.case.samples, named.case,
                    colname.control, number.control.samples, named.control, number.logFC, number.FDR) {
  
  # Load count data
  rawcounts <- read.csv(raw.reads.csv, header = TRUE, row.names = 1)
  
  # Identify case and control columns based on input parameters
  case_start <- which(colnames(rawcounts) == colname.case)
  case_end <- case_start + number.case.samples - 1
  control_start <- which(colnames(rawcounts) == colname.control)
  control_end <- control_start + number.control.samples - 1
  
  # Subset the count data to include only case and control columns
  counts <- rawcounts[, c(case_start:case_end, control_start:control_end)]
  
  # Create a DGEList object
  dge <- DGEList(counts = counts)
  
  # Filter lowly expressed genes
  keep <- rowSums(cpm(dge) > 0) >= ceiling(ncol(counts) / 2)
  dge <- dge[keep, ]
  
  # Perform normalization using TMM method
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Generate normalized data and save it
  norm_data <- cpm(dge, normalized.lib.sizes = TRUE)
  write.csv(norm_data, file = paste("normal_data_", named.case, "_", named.control, ".csv"), quote = TRUE, row.names = TRUE)
  
  # Create metadata for the samples
  metadata <- data.frame(
    sample_id = colnames(counts),
    condition = c(rep(named.case, times = number.case.samples), rep(named.control, times = number.control.samples))
  )
  metadata$condition <- relevel(factor(metadata$condition), ref = named.control)
  
  # Create design matrix
  design <- model.matrix(~ condition, data = metadata)
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit the model and perform likelihood ratio test for differential expression
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
  
  # Extract differentially expressed genes
  de_genes <- topTags(lrt, n = 99999999)$table
  de_genes_sig <- arrange(de_genes[de_genes$PValue < 0.05, ], logFC)
  
  # Add regulation information (up or down in the case group)
  tibble_data <- de_genes_sig %>% mutate(regulation = ifelse(logFC > 0, paste("up in", named.case), paste("down in", named.case)))
  
  # Convert to data frame and filter based on logFC and FDR thresholds
  result <- as.data.frame(tibble_data)
  filtered <- subset(result, abs(logFC) > number.logFC & FDR < number.FDR)
  
  # Save the filtered DEGs to a CSV file
  write.csv(filtered, file = paste(named.case, named.control, "logFC", number.logFC, "FDR", number.FDR, "DEGsNEW.csv", sep = "_"), quote = TRUE, row.names = TRUE)
}

# Function to set up argparse and pass arguments
run_with_args <- function() {
  parser <- ArgumentParser(description = "Differential Expression Analysis using edgeR")
  
  # Add arguments
  parser$add_argument("--raw.reads.csv", required = TRUE, help = "Path to the raw counts CSV file")
  parser$add_argument("--colname.case", required = TRUE, help = "Column name for the first case sample")
  parser$add_argument("--number.case.samples", type = "integer", required = TRUE, help = "Number of case samples")
  parser$add_argument("--named.case", required = TRUE, help = "Name for the case group")
  parser$add_argument("--colname.control", required = TRUE, help = "Column name for the first control sample")
  parser$add_argument("--number.control.samples", type = "integer", required = TRUE, help = "Number of control samples")
  parser$add_argument("--named.control", required = TRUE, help = "Name for the control group")
  parser$add_argument("--number.logFC", type = "numeric", required = TRUE, help = "logFC threshold")
  parser$add_argument("--number.FDR", type = "numeric", required = TRUE, help = "FDR threshold")
  
  # Parse the arguments
  args <- parser$parse_args()
  
  # Call the DEG function with parsed arguments
  get_DEG(args$raw.reads.csv, args$colname.case, args$number.case.samples, args$named.case,
          args$colname.control, args$number.control.samples, args$named.control,
          args$number.logFC, args$number.FDR)
}

# Run the script with arguments
run_with_args()

