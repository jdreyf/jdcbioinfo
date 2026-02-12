#' Plot histograms of significance (p-value & FDR) columns
#'
#' Plot histograms of significance (p-value & FDR) and QQ (p-value) columns.
#'
#' @param p.suffix Suffix for p-value columns. P-value column names cannot be duplicated.
#' @param fdr.suffix Suffix for FDR columns. Set to \code{NA} if no FDR columns. FDR column names cannot be duplicated.
#' @param sep Separator for column names before \code{p} or \code{FDR}, passed to \code{\link{extract_prefix}}.
#' If not found, it is assumed to be \code{NA}.
#' @param pi0 Logical indicating if proportion of null hypotheses should be calculated per p-value histogram.
#' @param method Method for estimating pi0, passed to \code{\link[limma:propTrueNull]{limma::propTrueNull}}. Ignored if \code{!pi0}.
#' @inheritParams ezlimmaplot::ezheat
#' @inheritParams ezlimmaplot::ezvenn
#' @details Some p-value columns must be identifiable using \code{p.suffix}. If \code{!is.na(fdr.suffix)}, FDR
#' colnames must have same prefix.
#' @return List of ggplots.
#' @export

# assume each comparison has a p-value & q-value column unless fdr.suffix=NA
# could allow for no prefix, ie colnames(tab)=c("p", "FDR")
signif_hist_qq <- function(tab, p.suffix="p", fdr.suffix="FDR", sep=".", pi0 = FALSE,
                           method=c("lfdr","convest"), name="signif_hist_qq"){
  stopifnot(nrow(tab) > 0, ncol(tab) > 0, !is.null(colnames(tab)))
  method <- match.arg(method)
  prefix.v <- extract_prefix(colnames(tab), suffix=p.suffix, sep=sep)

  if (any(duplicated(prefix.v))) stop("p-value column names are duplicated.")

  if (is.na(prefix.v[1])){
    p.cols <- match(p.suffix, colnames(tab))
  } else {
    p.cols <- match(paste0(prefix.v, sep, p.suffix), colnames(tab))
  }

  if (!is.na(fdr.suffix)){
    if (is.na(prefix.v[1])){
      fdr.cols <- match(fdr.suffix, colnames(tab))
    } else {
      fdr.cols <- match(paste0(prefix.v, sep, fdr.suffix), colnames(tab))
    }
    if (any(is.na(fdr.cols))) stop("!is.na(fdr.suffix) but FDR columns not found.")
    tab.ss <- tab[,sort(c(p.cols, fdr.cols))]
  } else {
    tab.ss <- tab[,p.cols]
  }

  plist <- list()
  for (ind.tmp in 1:length(p.cols)){
    prefix <- prefix.v[ind.tmp]
    p.col <- p.cols[ind.tmp]
    stopifnot(length(p.col)==1)
    subtitle <- NULL
    if(pi0){
      if (!requireNamespace("limma", quietly = TRUE)){
        stop("Package 'limma' needed to estimate pi0. Please install it.", call. = FALSE)
      }
      prop.null <- limma::propTrueNull(tab[,p.col], method = method)
      subtitle <- paste("Proportion of True Null = ", signif(prop.null, 3))
    }

    pvals <- tab[,p.col, drop=TRUE]
    ynull <- 0.05*length(pvals[!is.na(pvals)])
    p1 <- ggplot2::ggplot(tab, ggplot2::aes(x=!!rlang::sym(p.col))) +
      ggplot2::geom_histogram(bins=20) +
      ggplot2::labs(x="P-value", title=prefix, subtitle=subtitle) +
      ggplot2::geom_hline(yintercept=ynull, linetype="dashed") +
      ggplot2::annotate("text", x=0.8, y=ynull, label="Random", vjust=-1) +
      ggplot2::theme_minimal()


    tab <- tab %>%
      dplyr::mutate(expected=-log10((rank(!!rlang::sym(p.col), na.last="keep")-0.5)/sum(!is.na(!!rlang::sym(p.col)))),
                    observed=-log10(!!rlang::sym(p.col))) %>%
      dplyr::filter(!is.na(observed))
    lab <- expression("-"*log[10]~p*"-"*value)

    # Lambda GC (1 df)
    pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals>0 & pvals<=1]
    chisq_obs <- stats::qchisq(1-pvals, df=1)
    lambda_gc <- median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)

    p2 <- ggplot2::ggplot(tab.pv, ggplot2::aes(x=expected, y=observed)) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed") +
      ggplot2::labs(x=paste("Expected", lab),
                    y=paste("Observed", lab),
                    title=prefix) +
      ggplot2::annotate("text", x=0.2*max(tab.pv$expected), y=0.8*max(tab.pv$observed),
               label=sprintf("\u03BBGC = %.3f", lambda_gc), hjust=0) +
      ggplot2::theme_minimal()


    if (!is.na(fdr.suffix)){
      fdr.col <- fdr.cols[ind.tmp]
      stopifnot(length(fdr.col)==1)

      p3 <- ggplot2::ggplot(tab, ggplot2::aes(x=!!rlang::sym(fdr.col))) +
        ggplot2::geom_histogram(bins=20) +
        ggplot2::labs(x="FDR", title=prefix) +
        ggplot2::theme_minimal()
    }
    plist[[prefix]] <- patchwork::wrap_plots(p1, p2, if (!is.na(fdr.suffix)) p3 else NULL, ncol=2+!is.na(fdr.suffix))

  }

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"))
    on.exit(grDevices::dev.off())

    for (prefix in names(plist)){
      print(plist[[prefix]])
    }
  }

  return(invisible(plist))
}
