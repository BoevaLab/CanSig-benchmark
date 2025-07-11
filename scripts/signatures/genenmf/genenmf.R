


run_genenmf <- function(tumour_programs, n_cluster){
  library("logger")
  cc_genes <- c(Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes)
  cc_genes_threshold <- 10

    for (i in seq(n_cluster, n_cluster+10)){
      log_info("Trying ", i, " clusters.")
      results <- GeneNMF::getMetaPrograms(tumour_programs, nprograms = i, min.confidence = 0.0, max.genes = 50)
      metaprograms <- results$metaprograms.genes
      n_cc_genes <- list()
      for(name in names(metaprograms)){
          n_cc_genes[[name]] <- length(intersect(metaprograms[[name]], cc_genes))
      }
      if (sum(n_cc_genes < cc_genes_threshold) >= n_cluster){
        log_info("Found enough metasigs with fewer than ", cc_genes_threshold, " cell cycle genes.")
        break
      }
  }

metaprograms <- metaprograms[n_cc_genes < cc_genes_threshold]
metasigs <- data.frame(lapply(metaprograms, function(x) {
      x <- unlist(x)
      length(x) <- max(lengths(metaprograms))
      return(x)
    }))

colnames(metasigs) <- sprintf("Sig. %d", 1:ncol(metasigs))

    return(metasigs)
}


main <- function (){
  library("optparse")
  library("logger")
  library("scalop")
  option_list <- list(
          make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
          make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
          make_option(c("-n", "--n-cluster"), type="integer", default=NA, help="number of clusters"),
          make_option(c("-r", "--random-seed"), type="integer", default=NA)
          )

  opt <- parse_args(OptionParser(option_list=option_list))
  set.seed(opt$r)
  log_info(paste("Reading tumor programs from", opt$i))
  tumour_programs <- readRDS(opt$i)
  metasigs <- run_genenmf(tumour_programs = tumour_programs, n_cluster = opt$n)
  write.csv(metasigs, opt$o, row.names = F)
}

main()