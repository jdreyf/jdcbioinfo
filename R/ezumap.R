#' UMAP plot of first two dimensions
#'
#' UMAP plot of first two dimensions using \pkg{ggplot2}.
#'
#' @inheritParams ezlimmaplot::ezpca
#' @inheritParams ezlimmaplot::eztsne
#' @inheritParams umap::umap
#' @param ... list of settings with values overwrite defaults from UMAP \code{config} or passed to \code{\link[ggplot2:aes_]{aes_string}}.
#' @details \code{object} must have colnames, and if \code{pheno.df}
#' is given, it is checked that \code{colnames(object)==rownames(pheno.df)}.
#' @return Invisibly, a \code{ggplot} object. Its \code{data} element contains the first two principal components
#' appended to \code{pheno.df}.
#' @export

ezumap<- function(object, pheno.df, name='umap', pca=TRUE, initial_dims=nrow(pheno.df), config=umap::umap.defaults,
                  method=c("naive", "umap-learn"), preserve.seed=TRUE,
                  alpha=1, all.size=NULL, facet=NULL, title=NULL, subtitle=NULL, rm.leg.title=FALSE, labels=FALSE,
                  manual.color=NULL, manual.shape=NULL, plot=TRUE, ...){
  if (!requireNamespace("umap", quietly = TRUE)) stop("Package \"umap\" must be installed to use this function.", call. = FALSE)
  stopifnot(ncol(object)==nrow(pheno.df), colnames(object)==rownames(pheno.df))

  if(pca) {
    pca1 <- stats::prcomp(t(object[rowSums(is.na(object)) == 0, ]), scale. = FALSE)
    object <- pca1$x[, 1:initial_dims]
    }

  umap1 <- umap::umap(object, config=config, method = method, preserve.seed=preserve.seed, ...)
  umap1 <- as.data.frame(umap1$layout)
  colnames(umap1) <- c("UMAP1", "UMAP2")
  dat <- data.frame(pheno.df, umap1)

  dots <- list(...)
  if (is.null(names(dots))){
    n <- 0
  } else {
    chars <- vector("list", 2*length(dots))
    for(i in seq_along(dots)){
      chars[[2*i]] <- dots[[i]]
      chars[[2*i-1]] <- as.character(dat[, dots[[i]]])
    }
    n <- max(nchar(unlist(chars)))
  }

  width <- 7 + n / 12
  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"), width=width, height=7)
    on.exit(grDevices::dev.off())
  }

  #need to set alpha/all.size in geom_point, else it appears in legend
  qp <- ggplot2::ggplot(dat, mapping=ggplot2::aes_string(x='UMAP1', y='UMAP2', ...)) + ggplot2::theme_bw()

  if (!is.null(all.size)){
    qp <- qp + ggplot2::geom_point(size=all.size, alpha=alpha)
  } else {
    qp <- qp + ggplot2::geom_point(alpha=alpha)
  }

  if (!is.null(facet)){ qp <- qp + ggplot2::facet_grid(facet) }

  if (rm.leg.title){ qp <- qp + ggplot2::theme(legend.title=ggplot2::element_blank()) }

  if (!is.null(title)) { qp <- qp + ggplot2::ggtitle(label=title, subtitle=subtitle) }

  if (labels){
    dat2 <- dat
    dat2$row_names <- rownames(pheno.df)
    qp <- qp + ggplot2::geom_text(data=dat2, mapping=ggplot2::aes_string(label='row_names'), size=2, vjust=-.7)
  }

  if(!is.null(manual.color)) qp <- qp + ggplot2::scale_colour_manual(values=manual.color)
  if(!is.null(manual.shape)) qp <- qp + ggplot2::scale_shape_manual(values=manual.shape)

  if (plot) graphics::plot(qp)

  return(invisible(qp))
}
