library(multiMiR)
library(optparse)

# Define the function
miRNA_data_extractor <- function(input_file, species) {
  # Read the input file
  MM <- read.csv(input_file)
  
  # Extract unique gene IDs (previously miRNA names)
  unique_gene_ids <- unique(unlist(MM$geneid))
  
  # Loop over each gene ID and retrieve data
  for (i in unique_gene_ids) {
    tryCatch({
      # Fetch multiMiR data
      case1 <- get.multimir(org = species, target = i, table = "all", 
                            predicted.cutoff.type = "p", predicted.site = "all", 
                            summary = TRUE)
      data <- as.data.frame(case1@summary)
      
      # Write the output to a CSV file
      write.csv(data, paste(i, ".csv", sep = ""), row.names = FALSE)
      
    }, error = function(e) {
      # Print error message for failed miRNA retrieval
      print(paste("An error occurred for", i, ":", e))
    })
  }
}

# Command-line arguments setup
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input CSV file with gene IDs", metavar = "FILE"),
  make_option(c("-s", "--species"), type = "character", default = "mmu",
              help = "Species code for multiMiR (e.g., 'hsa' for human, 'mmu' for mouse)", 
              metavar = "SPECIES")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input file is provided
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("You must provide an input file with gene IDs.", call. = FALSE)
}

# Call the function with the provided arguments
miRNA_data_extractor(opt$input, opt$species)
