#' Apply `bacon` to z-scores
#'
#' Apply `bacon` function from package `bacon` to z-scores to obtain new z-scores, p-values, & FDR.
#'
#' @param tab Matrix-like object with statistical columns, some containing z-scores.
#' @param z.cols Indices or \code{\link{regexp}} with column names or column names suffix of numeric z-score columns.
#' @importFrom magrittr %>%
#' @return Data frame with bacon columns (z, p, FDR) inserted.
#' @export

# have not yet implemented alternative argument
ezbacon <- function(tab, z.cols="z"){
  if (!requireNamespace("bacon", quietly = TRUE)) stop("Package \"bacon\" must be installed to use this function.", call. = FALSE)
  stopifnot(ncol(tab) >= 1, nrow(tab) >= 1, !is.null(colnames(tab)))
  rnms_tab <- rownames(tab)
  z.colnms <- ezlimma:::grep_cols(tab, stat.cols=z.cols)
  # tab.z <- as.matrix(tab[, z.colnms])
  # if (!is.numeric(tab.z)) stop("Some z-score columns of tab are not numeric.")

  for (z.colnm in z.colnms){
    bc <- bacon::bacon(teststatistics = tab[, z.colnm, drop=FALSE])
    bacon.z <- bacon::tstat(bc)[,1]
    bacon.p <- 2*stats::pnorm(-abs(bacon.z))
    bacon.FDR <- stats::p.adjust(bacon.p, method="BH")

    prefix <- gsub(pattern = "\\..+", replacement = "", z.colnm)
    tab.bc <- data.frame(bacon.z, bacon.p, bacon.FDR)
    colnames(tab.bc) <- paste(prefix, colnames(tab.bc), sep="_")
    # strips row names
    tab <- tab %>% tibble::add_column(tab.bc, .after = z.colnm)
  }
  rownames(tab) <- rnms_tab
  tab
}
