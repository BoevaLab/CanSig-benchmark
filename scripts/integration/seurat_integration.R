library("Seurat")
library('optparse')


SAMPLE <- "sample"

run_integration <- function (data, type) {
  options(future.globals.maxSize = 8000 * 1024^2)
  k.weight <- min(100, min(plyr::count(data@meta.data[SAMPLE])["freq"]))

  batch_list = SplitObject(data, split.by = SAMPLE)
  features <- rownames(data)
  for (n_batch in 1:length(batch_list)){
    logger::log_info("Scaling data and computing PCA per sample.")
    batch_list[[n_batch]] <- ScaleData(batch_list[[n_batch]], features=features, verbose = FALSE)
    batch_list[[n_batch]] <- RunPCA(batch_list[[n_batch]], features=features, verbose = FALSE)
  }
  logger::log_info("Finding integration anchors.")
  anchors <- FindIntegrationAnchors(
    object.list = batch_list,
    anchor.features = rownames(data),
    reduction = type,
    l2.norm = TRUE,
    dims = 1:30,
    scale = FALSE,
    verbose = FALSE
  )
  logger::log_info("Integrating data.")

  integrated = IntegrateData(
              anchorset = anchors,
              new.assay.name = "integrated",
              k.weight=k.weight,
              verbose=FALSE
              )

  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated)
}


load_data <- function (data_path, n_hvg, batch_key){
  sce <- readRDS(data_path)
  features <- rownames(sce)[sce[["RNA"]]@meta.features[[paste0("highly_variable_", n_hvg, "_", batch_key)]]]
  sce <- subset(sce, features=features)
}

option_list <- list(
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
        make_option(c("-t", "--type"), type="character", default=NA),
        make_option(c("-r", "--random-seed"), type="integer", default=NA),
        make_option(c("--n-hvg"), type="character", default=NA),
        make_option(c("--batch-key"), type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$r)
data <- load_data(data_path = opt$i, n_hvg=opt$n_hvg, batch_key=opt$batch_key)
data <- subset(data, subset = phase=="G1")
integrated <- run_integration(data, type=opt$t)

if (grepl("\\.csv$", opt$o, ignore.case = TRUE)) {
  logger::log_info("Saving integrated latent representation as CSV.")
  latent_repr <- integrated@reductions$pca@cell.embeddings
  latent_df <- as.data.frame(latent_repr)
  write.csv(latent_df, file = opt$o, row.names = TRUE)
} else {
  logger::log_info("Saving integrated object as RDS.")
  saveRDS(integrated, file = opt$o)
}

warnings()
