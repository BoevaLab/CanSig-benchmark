library('optparse')
library('sceasy')

option_list <- list(
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file")
)


opt = parse_args(OptionParser(option_list=option_list))


.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[["name"]] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        paste("Dropping single category variables:"),
        paste(colnames(df)[k_singular], collapse = ", ")
      )
    }
    df <- df[, !k_singular, drop = F]
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
  }
  return(df)
}


.obs2metadata <- function(obs_pd, assay = "RNA") {
  obs_df <- .regularise_df(reticulate::py_to_r(obs_pd), drop_single_values = FALSE)
  attr(obs_df, "pandas.index") <- NULL
  return(obs_df)
}


.var2feature_metadata <- function(var_pd) {
  var_df <- .regularise_df(reticulate::py_to_r(var_pd), drop_single_values = FALSE)
  attr(var_df, "pandas.index") <- NULL
  return(var_df)
}


.layer_to_matrix <- function (layer, colnames, rownames){
  sp <- reticulate::import("scipy.sparse", convert = FALSE)
  if (reticulate::py_to_r(sp$issparse(layer))) {
    X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(layer)))
  } else {
    X <- t(reticulate::py_to_r(layer))
  }
  X <- as(X, "dgCMatrix")
  colnames(X) <- colnames
  rownames(X) <- rownames
  return(X)
}



anndata2seurat <- function(opt, assay = "RNA") {
  COUNT_LAYER <- "counts"
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  inFile <- path.expand(opt$i)

  anndata <- reticulate::import("anndata", convert = FALSE)

  ad <- anndata$read_h5ad(inFile)

  obs_df <- .obs2metadata(ad$obs)
  var_df <- .var2feature_metadata(ad$var)

  log_tpm <- .layer_to_matrix(ad$X, rownames(obs_df), rownames(var_df))
  X <- .layer_to_matrix(ad$layers[[COUNT_LAYER]], rownames(obs_df), rownames(var_df))

  assays <- list(Seurat::CreateAssayObject(counts = X))
  assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = "data", new.data = log_tpm)
  message("X -> raw.data; layers[<tpm_key>] -> data")

  names(assays) <- assay
  Seurat::Key(assays[[assay]]) <- paste0(tolower(assay), "_")
  assays[[assay]]@meta.features <- var_df


  project_name <- sub("\\.h5ad$", "", basename(inFile))
  srt <- new("Seurat", assays = assays, project.name = project_name, version = packageVersion("Seurat"))
  Seurat::DefaultAssay(srt) <- assay
  Seurat::Idents(srt) <- project_name

  srt@meta.data <- obs_df
  embed_names <- unlist(reticulate::py_to_r(ad$obsm_keys()))
  if (length(embed_names) > 0) {
    embeddings <- sapply(embed_names, function(x) as.matrix(reticulate::py_to_r(ad$obsm[x])), simplify = FALSE, USE.NAMES = TRUE)
    names(embeddings) <- embed_names
    for (name in embed_names) {
      rownames(embeddings[[name]]) <- colnames(assays[[assay]])
    }

    dim.reducs <- vector(mode = "list", length = length(embeddings))
    for (i in seq(length(embeddings))) {
      name <- embed_names[i]
      embed <- embeddings[[name]]
      key <- switch(name,
        sub("_(.*)", "\\L\\1", sub("^X_", "", toupper(name)), perl = T),
        "X_pca" = "PC",
        "X_tsne" = "tSNE",
        "X_umap" = "UMAP"
      )
      colnames(embed) <- paste0(key, "_", seq(ncol(embed)))
      dim.reducs[[i]] <- Seurat::CreateDimReducObject(
        embeddings = embed,
        loadings = new("matrix"),
        assay = assay,
        stdev = numeric(0L),
        key = paste0(key, "_")
      )
    }
    names(dim.reducs) <- sub("X_", "", embed_names)

    for (name in names(dim.reducs)) {
      srt[[name]] <- dim.reducs[[name]]
    }
  }


  saveRDS(object = srt, file = opt$o)
}


anndata2seurat(opt)