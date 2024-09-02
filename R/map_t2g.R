#' Map pathway term to gene IDs to another ID type
#'
#' Map pathway term to gene with old gene IDs (e.g. entrez IDs) to another ID type (e.g. gene symbols) using a table with both
#' annotation types. The original IDs can be separated within a column of annot with a string, e.g. " /// ".
#' The other IDs should be the rownames of annot. The idea for this is based on the Gene Set Enrichment
#' Analysis (GSEA) utility [chip2chip](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?Chip2Chip_Page).
#'
#' @param t2g TERM TO GENE mapping, a data.frame of 2 column with term and gene
#' @param annot annotation for the features that has a column of the same type as in gene set list *G*.
#' @param sep.str strings that separates symbols if there are multiple symbols for a feature.
#' @param symbol.col column name or index for the symbol column in *annot*.
#' @param fixed logical passed to \code{strsplit}; if \code{TRUE}, match \code{sep.str} exactly, otherwise use regular expression.
#' @inheritParams roast_contrasts
#' @return TERM TO GENE data.frame with the symbols replaced by the mapped IDs.
#' @details Gene annotations are transformed to upper-case to avoid missing matches of different cases.
#' @export

map_t2g <- function(t2g, annot, sep.str=" /// ", symbol.col="Gene.Symbol", fixed=FALSE){

  stopifnot(symbol.col %in% colnames(annot) || symbol.col %in% 1:ncol(annot), !is.null(rownames(annot)))

  sym.v <- toupper(annot[, symbol.col, drop = TRUE])
  sym.lst <- strsplit(sym.v, split = sep.str, fixed=fixed)
  row.nms <- rep(rownames(annot), times = vapply(sym.lst, FUN=length, FUN.VALUE = numeric(1)))
  map <- data.frame(Rowname = row.nms, Symbol = unlist(sym.lst), stringsAsFactors = FALSE)

  t2g <- dplyr::group_split(t2g, !!rlang::sym(colnames(t2g)[1]))
  t2g <- lapply(t2g, FUN = as.data.frame)

  mapped.t2g <- list()
  for(i in seq_along(t2g)){
    gene <- unique(map$Rowname[ map$Symbol %in% toupper(t2g[[i]][,2]) ])

    if (length(gene) == 0) {
      mapped.t2g[[i]] <- NULL
    } else {
      mapped.t2g[[i]] <- data.frame(term = unique(t2g[[i]][,1]), gene = gene)
    }
  }
  mapped.t2g <- Reduce(rbind, mapped.t2g)

  return(mapped.t2g)
}
