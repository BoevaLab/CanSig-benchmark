library("Seurat")
library('optparse')
library('logger')

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


runFastMNN = function(sobj, features) {
  suppressPackageStartupMessages({
    require(batchelor)
  })

  if (is.null(sobj@assays$RNA)) {
    # Seurat v4
    expr <- GetAssayData(sobj, slot = "data")
  } else {
    # Seurat v3
    expr <- sobj@assays$RNA@data
  }

  sce <- fastMNN(expr, batch = sobj@meta.data[["sample"]], subset.row=features)
  corrected_data <- as.matrix(assay(sce, "reconstructed"))
  
  sobj[["integrated"]] <- CreateAssayObject(data=corrected_data)
  DefaultAssay(sobj) <- "integrated"
  logger::log_info("Scaling integrated data")
  sobj <- ScaleData(sobj, verbose = FALSE)
  logger::log_info("Running PCA on integrated data")
  sobj <- RunPCA(sobj, verbose = FALSE, features=rownames(sobj))
  return(sobj)
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
  integrated <- runFastMNN(data, features = features)
  
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



main()