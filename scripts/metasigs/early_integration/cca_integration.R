

library("Seurat")
library('optparse')
SAMPLE <- "sample"

run_integration <- function (data, k.weight, hvg=4000) {

  batch_list = SplitObject(data, split.by = SAMPLE)
  for (n_batch in 1:length(batch_list)){
    batch_list[[n_batch]] <- ScaleData(batch_list[[n_batch]])
  }

  anchors = FindIntegrationAnchors(
              object.list = batch_list,
              anchor.features = hvg,
              l2.norm = T,
              scale = T,
              dims = 1:30,
              k.anchor = 5,
              k.filter = 200,
              k.score = 30,
              max.features = 200,
              eps = 0)
  integrated = IntegrateData(
              anchorset = anchors,
              new.assay.name = "integrated",
              k.weight = k.weight)


  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated)
}


load_data <- function (data_path){
  data <- readRDS(data_path)
}

option_list <- list(
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file")
)

opt <- parse_args(OptionParser(option_list=option_list))
data <- load_data(data_path = opt$i)
data <- subset(data, subset = phase=="G1")

k.weight = min(100, min(plyr::count(data@meta.data[SAMPLE])["freq"]))
integrated <- run_integration(data, k.weight=k.weight)
saveRDS(integrated, file = opt$o)



