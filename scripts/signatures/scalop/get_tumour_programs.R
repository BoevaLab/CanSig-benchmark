get_highly_expressed_genes <- function(sce){

  data <- Seurat::GetAssayData(sce, "data")
  data@x <- 2^(data@x) - 1
  gene_mean_expression <- log2(Matrix::rowMeans(10*data)+1)
  return(names(gene_mean_expression)[gene_mean_expression>4.])
}

 get_programs <- function(m, groups = NULL, nsig1 = 50, nsig2 = 10, jaccard = 0.7, p = 0.01, lfc = log2(2), pmethod = "BH") {
  library(scalop)
  if (is.null(groups)) {
    groups = hca_groups(rowcenter(m))
  }
  deas = dea(m, group = groups, return.val = "df", p = p, lfc = lfc, arrange.by = "lfc", pmethod = pmethod)
  sig1 = sapply(deas, function(df) sum(df$p.adj <= 0.01, na.rm = TRUE))
  sig2 = sapply(deas, function(df) sum(df$p.adj <= 0.001, na.rm = TRUE))
  sig3 = sapply(deas, function(df) sum(df$p.adj <= 0.0001, na.rm = TRUE))
  bool = sig1 >= nsig1 & sig2 >= nsig2
  if (sum(bool) == 0) {
    return(NULL)
  }
  deas = deas[bool]
  groups = groups[bool]
  sig1 = sig1[bool]
  sig2 = sig2[bool]
  sig3 = sig3[bool]
  ord = order(sig1, sig2, sig3, decreasing = T)
  deas = deas[ord]
  groups = groups[ord]
  sig1 = sig1[ord]
  sig2 = sig2[ord]
  sig3 = sig3[ord]
  jac.pass = jacFilt(groups, threshold = jaccard, which = TRUE)
  deas = deas[jac.pass]
  groups = groups[jac.pass]
  sig1 = sig1[jac.pass]
  sig2 = sig2[jac.pass]
  sig3 = sig3[jac.pass]
  programs = sapply(deas, function(d) d$gene[1 : nsig1], simplify = F)
  lfcs = dea(m, group = groups, return.val = "lfc", arrange.by = "none", p = NULL, lfc = NULL, pmethod = pmethod)
  return(list(programs = programs, profiles = lfcs, groups = groups, deas = deas, sig1 = sig1, sig2 = sig2, sig3 = sig3))
}

get_tumor_programs <- function(sce_list, row_means){
  prog.obj <- list()
  for(batch_key in names(sce_list)){
    log_info(paste("Processing batch:", batch_key))
    batch <- sce_list[[batch_key]]
    input <- t(scale(t(as.matrix(Seurat::GetAssayData(batch, "data"))), center = row_means, scale = F))
    groups <- scalop::hca_groups(input, max.size = 0.5)
    res <- get_programs(input, groups=groups , jaccard = 0.75, p=0.05, lfc = log2(3))
    if(is.null(res)){
      log_info(paste0("No differentially expresed groups found for: ", batch_key))
      next
    }
    prog.obj <- c(prog.obj, setNames(list(res), batch_key))
  }


  names(prog.obj) <- paste0(names(prog.obj), "..")

  tumour_programs <- sapply(prog.obj, `[[`, 'programs', simplify = F) %>%
    unlist(., recursive = F)
  return(tumour_programs)
}


get_row_means <- function(sce){
  return(Matrix::rowMeans(Seurat::GetAssayData(sce, "data")))
}


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
  log_info(paste("Reading data from", opt$i))
  sce <- readRDS(opt$i)

  selected_genes <- get_highly_expressed_genes(sce)
  log_info(paste("Running with", length(selected_genes), "genes"))
  sce <- subset(sce, features = selected_genes)
  row_means <- get_row_means(sce)
  sce <- Seurat::SplitObject(sce, "sample")

  tumour_programs <- get_tumor_programs(sce_list = sce, row_means=row_means)
  log_info("Found", length(tumour_programs), "programs.")
  log_info("Saving rds to", opt$o)
  saveRDS(tumour_programs, opt$o)
}

main()
