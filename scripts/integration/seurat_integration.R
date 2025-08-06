library("Seurat")
library('optparse')
library('logger')

# Constants
SAMPLE_COLUMN <- "sample"
INTEGRATION_DIMS <- 1:30
MAX_FUTURE_GLOBALS_SIZE <- 1000000 * 1024^2

#' Run integration on Seurat data
#' @param data Seurat object
#' @param features Character vector of feature names for integration
#' @param reduction_type Type of reduction to use for integration
#' @return Integrated Seurat object
run_integration <- function(data, features, reduction_type) {
  logger::log_info("Starting data integration process")
  
  # Validate inputs
  if (is.null(data) || !inherits(data, "Seurat")) {
    logger::log_error("Invalid Seurat object provided")
    stop("Data must be a valid Seurat object")
  }
  
  if (length(features) == 0) {
    logger::log_error("No features provided for integration")
    stop("Features list cannot be empty")
  }
  
  logger::log_info("Setting future globals max size to {MAX_FUTURE_GLOBALS_SIZE} bytes")
  plan("multicore", workers = 5)
  options(future.globals.maxSize = MAX_FUTURE_GLOBALS_SIZE)
  
  logger::log_info("Splitting data by sample column: {SAMPLE_COLUMN}")
  batch_list <- SplitObject(data, split.by = SAMPLE_COLUMN)
  logger::log_info("Created {length(batch_list)} batches for integration")
  
  for (n_batch in 1:length(batch_list)){
  logger::log_info("Scaling data and running PCA per sample.")
    batch_list[[n_batch]] <- ScaleData(batch_list[[n_batch]], features=features, verbose = FALSE)
    batch_list[[n_batch]] <-  RunPCA(batch_list[[n_batch]], features=features, verbose = FALSE)
  }

  logger::log_info("Finding integration anchors using {length(features)} features")
  logger::log_debug("Integration parameters: reduction={reduction_type}, dims={paste(INTEGRATION_DIMS, collapse=',')}")

  anchors <- FindIntegrationAnchors(
    object.list = batch_list,
    anchor.features = features,
    reduction = reduction_type,
    l2.norm = TRUE,
    dims = INTEGRATION_DIMS,
    normalization.method = "LogNormalize",
    scale = TRUE,
    verbose = FALSE
  )
  
  # Calculate k.weight based on anchor distribution
  k_weight <- min(100, min(table(anchors@anchors[5])))
  logger::log_info("Calculated k.weight: {k_weight}")

  logger::log_info("Integrating data using anchors")
  integrated <- tryCatch({
    IntegrateData(
      anchorset = anchors,
      new.assay.name = "integrated",
      k.weight = k_weight,
      verbose = FALSE
    )
  }, error = function(e) {
    logger::log_warn("Integration failed with k.weight={k_weight}: {e$message}")
    logger::log_info("Retrying integration with k.weight=30")
    
    tryCatch({
      IntegrateData(
        anchorset = anchors,
        new.assay.name = "integrated",
        k.weight = 30,
        verbose = FALSE
      )
    }, error = function(e2) {
      logger::log_error("Integration failed even with k.weight=30: {e2$message}")
      stop("Data integration failed with both calculated and fallback k.weight values")
    })
  })
  
  logger::log_info("Scaling integrated data")
  #integrated <- ScaleData(integrated, verbose = FALSE)
  
  logger::log_info("Running PCA on integrated data")
  integrated <- RunPCA(integrated, verbose = FALSE)
  
  logger::log_info("Integration process completed successfully")
  return(integrated)
}

#' Load Seurat data from RDS file
#' @param data_path Path to RDS file containing Seurat object
#' @return Seurat object
load_data <- function(data_path) {
  logger::log_info("Loading data from: {data_path}")
  
  if (!file.exists(data_path)) {
    logger::log_error("Input file does not exist: {data_path}")
    stop("Input file not found")
  }
  
  if (!grepl("\\.rds$", data_path, ignore.case = TRUE)) {
    logger::log_warn("Input file does not have .rds extension, proceeding anyway")
  }
  
  tryCatch({
    data <- readRDS(data_path)
    logger::log_info("Successfully loaded data with {ncol(data)} cells and {nrow(data)} features")
    return(data)
  }, error = function(e) {
    logger::log_error("Failed to load data: {e$message}")
    stop("Could not read RDS file")
  })
}


# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), 
              type = "character", 
              default = NA, 
              help = "Input RDS file containing Seurat object"),
  make_option(c("-o", "--output"), 
              type = "character", 
              default = NA, 
              help = "Output file (.rds for Seurat object, .csv for latent representation)"),
  make_option(c("-t", "--type"), 
              type = "character", 
              default = NA,
              help = "Integration reduction type"),
  make_option(c("-r", "--random-seed"), 
              type = "integer", 
              default = NA,
              help = "Random seed for reproducibility"),
  make_option(c("--n-hvg"), 
              type = "character", 
              default = NA,
              help = "Number of highly variable genes parameter",
              dest = "n_hvg"),
  make_option(c("--batch-key"), 
              type = "character", 
              default = NA,
              help = "Batch key for highly variable genes selection",
              dest = "batch_key")
)

# Main execution
main <- function() {
  logger::log_info("=== Starting Seurat Integration Pipeline ===")
  
  # Parse command line arguments
  opt <- parse_args(OptionParser(option_list = option_list))
  logger::log_info("Arguments: input={opt$i}, output={opt$o}, type={opt$t}, seed={opt$r}, n_hvg={opt$n_hvg}, batch_key={opt$batch_key}")
  
  # Set random seed for reproducibility
  logger::log_info("Setting random seed to: {opt$r}")
  set.seed(opt$r)
  
  # Load data
  data <- load_data(data_path = opt$i)
  
  # Extract parameters
  n_hvg <- opt$n_hvg
  batch_key <- opt$batch_key
  
  # Select highly variable features
  hvg_column <- paste0("highly_variable_", n_hvg, "_", batch_key)
  logger::log_info("Selecting highly variable genes using column: {hvg_column}")
  
  if (!hvg_column %in% colnames(data[["RNA"]]@meta.features)) {
    logger::log_error("Highly variable genes column not found: {hvg_column}")
    stop("Invalid HVG column specification")
  }
  
  features <- rownames(data)[data[["RNA"]]@meta.features[[hvg_column]]]
  logger::log_info("Selected {length(features)} highly variable features")
  
  if (length(features) == 0) {
    logger::log_error("No highly variable features found")
    stop("No features selected for integration")
  }
  
  # Subset data to G1 phase cells
  logger::log_info("Subsetting data to G1 phase cells")
  original_cell_count <- ncol(data)
  
  if (!"phase" %in% colnames(data@meta.data)) {
    logger::log_error("Phase column not found in metadata")
    stop("Phase information required for subsetting")
  }
  
  data <- subset(data, subset = phase == "G1")
  final_cell_count <- ncol(data)
  
  logger::log_info("Filtered from {original_cell_count} to {final_cell_count} G1 cells")
  
  if (final_cell_count == 0) {
    logger::log_error("No G1 phase cells found after filtering")
    stop("No cells remaining after G1 filtering")
  }
  
  # Run integration
  integrated <- run_integration(data, features = features, reduction_type = opt$t)
  
  # Save results
  if (grepl("\\.csv$", opt$o, ignore.case = TRUE)) {
    logger::log_info("Saving integrated latent representation as CSV: {opt$o}")
    
    if (!"pca" %in% names(integrated@reductions)) {
      logger::log_error("PCA reduction not found in integrated object")
      stop("PCA reduction missing")
    }
    
    latent_repr <- integrated@reductions$pca@cell.embeddings
    latent_df <- as.data.frame(latent_repr)
    
    write.csv(latent_df, file = opt$o, row.names = TRUE)
    logger::log_info("Successfully saved {nrow(latent_df)} cell embeddings with {ncol(latent_df)} components")
    
  } else {
    logger::log_info("Saving integrated Seurat object as RDS: {opt$o}")
    saveRDS(integrated, file = opt$o)
    logger::log_info("Successfully saved integrated Seurat object")
  }
  
  logger::log_info("=== Integration pipeline completed successfully ===")
}

# Execute main function
if (!interactive()) {
  main()
}