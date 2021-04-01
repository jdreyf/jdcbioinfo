#' Map gene IDs in geneset list to rownames of annotation data frame
#'
#' Map gene IDs in geneset list to rownames of annotation data frame. The name is based on the Gene Set Enrichment
#' Analysis (GSEA) utility [chip2chip](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?Chip2Chip_Page).
#'
#' @param G Geneset list such as those read in with `ezlimma::read.gmt` with list elements whose 3rd element
#' are the gene IDs.
#' @param annot Annotation matrix or data frame with row names and column names that has a column with gene IDs t
#' hat match those in `G` and whose row names represent the new gene IDs we want.
#' @param geneid.col Column name or column index of `annot` with gene IDs that match those in `G`.
#' @param split.str Character vector (or object which can be coerced to such) containing regular expression
#' to use for splitting `annot[, geneid.col]` when it contains multiple separated gene IDs.
#' @details Both `rownames(annot)` and gene IDs in `G` are converted to upper case to aid matching.
#' @return `G` with new gene IDs matching `rownames(annot)`.
#' @export

ezchip2chip <- function(G, annot, geneid.col='Gene.Symbol', split.str=' /// '){
  stopifnot(geneid.col %in% colnames(annot) || geneid.col %in% 1:ncol(annot))
  sym.v <- toupper(annot[, geneid.col, drop=TRUE])
  sym.lst <- strsplit(sym.v, split = split.str)
  row.ind <- unlist(sapply(seq_along(sym.lst), FUN = function(i) rep(i, length(sym.lst[[i]])) ))
  id.map <- data.frame(Rowname = rownames(annot)[row.ind], Symbol = unlist(sym.lst))
  for(i in seq_along(G)){
    G[[i]][[3]] <- unique(id.map$Rowname[ id.map$Symbol %in% toupper(G[[i]][[3]]) ])
  }
  return(G)
}

