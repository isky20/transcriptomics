# Load necessary libraries
library(readxl)
library(multiMiR)

# Define the function
get_mirna_targets <- function(input_file, output_folder) {
  # Import the file
  MM <- read.csv(input_file)
  mirna <- MM$miRNA.names

  # Set working directory to the output folder
  setwd(output_folder)

  # Initialize a list to save miRNA targets
  target <- list()

  # Get miRNA targets using tryCatch to handle errors
  for(i in mirna) {
    tryCatch({
      target[[i]] <- get_multimir(mirna = noquote(i), summary = TRUE)
    }, error = function(e) {
      message(paste("Error for miRNA:", i, "-", e$message))
    })
  }

  # Save each miRNA's targets to a CSV file
  for(k in seq_along(target)) {
    print(names(target)[k])
    write.csv(target[[k]]@data, paste0(names(target)[k], ".csv"))
  }
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript get_mirna_targets.R <input_file> <output_folder>")
}

input_file <- args[1]
output_folder <- args[2]

# Call the function with specified parameters
get_mirna_targets(input_file, output_folder)
