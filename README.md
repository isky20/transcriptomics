### Load Required Libraries:

Loads necessary R packages (edgeR for differential expression analysis, dplyr for data manipulation).
### Define get_DEG Function:

Defines a function get_DEG that performs differential expression analysis using edgeR.
Parameters include file paths, column names, sample sizes, log fold change (logFC), and false discovery rate (FDR) thresholds.
### Load Count Data:

Reads raw count data from a CSV file (raw.reads.csv) into R using read.csv.
### Subset Data:

Identifies and subsets columns corresponding to case and control samples based on user-provided column names and sample sizes.
### Create DGEList:

Constructs a DGEList object (dge) from the subsetted count data for further analysis.
### Filter Lowly Expressed Genes:

Filters out genes with low expression levels using a count-per-million (CPM) threshold to ensure robust differential expression analysis.
### Normalization:

Normalizes the count data using the TMM (trimmed mean of M-values) method to adjust for differences in library sizes across samples.
### Save Normalized Data:

Writes the normalized count data (norm_data) to a CSV file for future reference.
### Create Metadata:

Generates metadata including sample IDs and conditions (case or control) for sample classification and statistical modeling.
### Create Design Matrix:

Constructs a design matrix (design) specifying the experimental conditions (case vs. control) for differential expression analysis.
### Estimate Dispersion:

Estimates dispersion of gene expression across samples using the fitted design matrix to prepare for statistical modeling.
### Differential Expression Analysis:

Fits a generalized linear model (glmFit) using the edgeR package to model differential gene expression between case and control samples.
Performs a likelihood ratio test (glmLRT) to assess the significance of differential expression.
### Extract DE Genes:

Extracts genes that are differentially expressed based on a specified significance threshold (adjusted p-value < 0.05).
### Add Regulation Information:

Adds information on gene regulation (up-regulated or down-regulated in the case group compared to the control group) based on log fold changes (logFC).
### Filter Based on LogFC and FDR:

Filters the list of differentially expressed genes (DEGs) based on user-defined thresholds for log fold change (number.logFC) and false discovery rate (number.FDR).
### Save Filtered DEGs:

Writes the filtered list of differentially expressed genes (filtered) to a CSV file with a descriptive name reflecting the conditions and thresholds used (named.case, named.control, number.logFC, number.FDR).
### Example Usage:

Demonstrates an example usage of the get_DEG function with sample parameters to perform differential expression analysis on RNA-seq data.
