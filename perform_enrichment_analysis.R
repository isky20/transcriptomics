#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readxl)
library(rbioapi)

# Define the function
perform_enrichment_analysis <- function(file_name, organism_id) {
  
  # Check if the file exists
  if (!file.exists(file_name)) {
    stop("File does not exist: ", file_name)
  }
  
  # Read CSV file into a data frame
  df <- read.csv(file_name)
  
  # Analyze up-regulated genes
  up_list <- df$Symbol[grepl("up", df$regulation, ignore.case = TRUE)]
  enrichmentup <- rbioapi::rba_string_enrichment(up_list, organism_id, split_df = FALSE)
  
  if (nrow(enrichmentup) == 0) {
    cat("No enrichment found for up-regulated genes in", file_name, "\n")
  } else {
    subset_dfUP <- enrichmentup[, c("term", "description", "category", "number_of_genes")]
    write.csv(subset_dfUP, paste0(gsub(".csv", "", file_name), "_up.csv"), row.names = FALSE)
  }
  
  # Analyze down-regulated genes
  down_list <- df$Symbol[grepl("down", df$regulation, ignore.case = TRUE)]
  enrichment <- rbioapi::rba_string_enrichment(down_list, organism_id, split_df = FALSE)
  
  if (nrow(enrichment) == 0) {
    cat("No enrichment found for down-regulated genes in", file_name, "\n")
  } else {
    subset_df <- enrichment[, c("term", "description", "category", "number_of_genes")]
    write.csv(subset_df, paste0(gsub(".csv", "", file_name), "_down.csv"), row.names = FALSE)
  }
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript perform_enrichment_analysis.R <file_name> <organism_id>")
}
file_name <- args[1]
organism_id <- as.numeric(args[2])

# Run the function with command-line arguments
perform_enrichment_analysis(file_name, organism_id)
