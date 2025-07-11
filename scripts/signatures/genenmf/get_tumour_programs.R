
main <- function (){
  library("optparse")
  library("logger")
  option_list <- list(
          make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
          make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
          make_option(c("-r", "--random-seed"), type="integer", default=NA)

          )

  opt <- parse_args(OptionParser(option_list=option_list))
  set.seed(opt$r)
  sce <- readRDS(opt$i)

  splitobj.list <- Seurat::SplitObject(sce, split.by="sample")
  tumour_programs <- GeneNMF::multiNMF(splitobj.list, k = 3:12, max.exp = 100000, nfeatures = 8000)
  log_info("Found", length(tumour_programs), "programs.")
  log_info("Saving rds to", opt$o)
  saveRDS(tumour_programs, opt$o)
}

main()
