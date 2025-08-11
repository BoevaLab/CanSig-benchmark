
library("Seurat")
library('optparse')
library('logger')
SAMPLE <- "sample"



binary_search_leiden <- function (data, n_cluster, start, end, random_seed, epsilon){
      while (end-start > epsilon){
        mid = (end + start) / 2
        data <- Seurat::FindClusters(data, resolution = mid, n.start=1, verbose = F, random.seed=random_seed)
        n_tmp = length(unique(data@active.ident))
        log_info("Found ", n_tmp, " cluster for resolution ", mid)
        if (n_tmp == n_cluster){
          return (data)
        }
        if (n_tmp > n_cluster){
            end = mid
        }
        if (n_tmp < n_cluster){
            start = mid
          }
    }
}

get_clusters <- function (data, n_cluster, random_seed,start = 1e-4, end = 2., epsilon = 1e-8){
  for (i in 1:10){
    log_info("Trying randomseed ", i, ".")
    results <- binary_search_leiden(data = data, n_cluster=n_cluster, start = start, end = end, epsilon = epsilon, random_seed = 10000*i+random_seed)
    if (!is.null(results)){
      return(results)
    }
  }
}

load_data <- function (data_path){
  data <- readRDS(data_path)
}

option_list <- list(
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
		    make_option(c("-n", "--n-cluster"), type="integer", default=NA, help="output file"),
        make_option(c("-r", "--random-seed"), type="integer", default=NA)


)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$r)
integrated <- load_data(opt$i)
integrated <-FindNeighbors(integrated, reduction = "pca")
integrated <-get_clusters(integrated, opt$n, opt$r)
DefaultAssay(integrated) <- "RNA"
results = data.frame(row.names = 1:50)
for (group in unique(integrated@active.ident)){
  if (sum(integrated@active.ident == group)< 5)
    {next}
  marker <- FindMarkers(integrated, ident.1 = group, only.pos = T)
  marker <- row.names(marker)[1:50]
  results[paste("Sig.", as.character(as.integer(group)+1))] <- marker
}

write.csv(results, opt$o, row.names = F)
