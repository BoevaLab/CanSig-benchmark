
get_programs <- function (tumour_programs, n_cluster){
    cc.matrix <- Jaccard(tumour_programs)

    cluster <- scalop::hca_groups(cc.matrix,
                                         cor.method="none",
                                         k=n_cluster,
                                         min.size=0,
                                         max.size=1)

    mp_freqs <- sapply(cluster, function(k) sort(table(unlist(tumour_programs[k])), decreasing = T), simplify = F)
    metaprograms <- sapply(mp_freqs, function(tab) head(names(tab), 50), simplify = F)
    return(list(metaprograms=metaprograms, cluster=cluster))
}


run_scalop <- function(tumour_programs, n_cluster){
  results <- get_programs(tumour_programs, 2)
    metaprograms <- results$metaprograms
    cluster_cellcycle <- results$cluster
    cc_genes <- c(Seurat::cc.genes$s.genes, Seurat::cc.genes$g2m.genes)

    log_info("Removing high cc signatures")

    n_cc_genes <- list()
    for(name in names(metaprograms)){
        n_cc_genes[[name]] <- length(intersect(metaprograms[[name]], cc_genes))
    }

    if(any(n_cc_genes > 0.25 * lengths(metaprograms))){
        noncc_cluster_names <- names(which.min(n_cc_genes))
        log_info("Found a meta-signature cluster associated with cell cylce")
        noncc_program_names <- cluster_cellcycle[[noncc_cluster_names]]
        log_info("Keeping ", length(noncc_program_names), " signatures.")
    }else{
        log_info("No meta-signature associated with cell cycle found. Retaining all sigantures.")
        noncc_program_names <- names(tumour_programs)
    }


  noncc_programs <- tumour_programs[noncc_program_names]
  results <- get_programs(noncc_programs, n_cluster = n_cluster)
  metaprograms <- results$metaprograms

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
  log_info("Loaded ", length(tumour_programs), "signatures.")

  metasigs <- run_scalop(tumour_programs = tumour_programs, n_cluster = opt$n)
  write.csv(metasigs, opt$o, row.names = F)
}

main()