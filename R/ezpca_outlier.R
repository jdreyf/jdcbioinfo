#' PCA plot with optional outlier detection
#'
#' PCA scatter plot of the first two principal components. Optionally draws a
#' confidence ellipse based on Mahalanobis distance and labels points outside
#' this ellipse as potential outliers. All original functionality of
#' \code{\link{ezlimmaplot::ezpca}} is preserved.
#'
#' @param outlier.pval Numeric p‑value threshold for outlier detection
#'   (e.g., 0.05). If \code{NULL} (default), no outlier ellipse is drawn and no
#'   outlier labels are added.
#' @param label.outliers Logical; if \code{TRUE} (default) and \code{outlier.pval}
#'   is not \code{NULL}, labels points outside the ellipse using \code{ggrepel}.
#'   The global \code{labels} parameter is automatically set to \code{FALSE}
#'   when outlier detection is active.
#' @inheritParams ezlimmaplot::ezpca
#' @param ... Aesthetics passed to \code{\link[ggplot2:aes_]{aes_string}}.
#'
#' @details PCA is computed with \code{\link[stats]{prcomp}} (no scaling).
#' The squared Mahalanobis distance for a point with PC scores \eqn{t_1, t_2} is
#' \eqn{D^2 = t_1^2/\lambda_1 + t_2^2/\lambda_2}, where \eqn{\lambda_1,\lambda_2}
#' are the variances (eigenvalues) of the first two components.
#' Under bivariate normality, \eqn{D^2 \sim \chi^2_2}. Points with
#' \eqn{D^2 > \chi^2_{2,1-p}} are considered outliers.
#'
#' When \code{outlier.pval} is provided, a single confidence ellipse (ignoring
#' any grouping aesthetics) is drawn using \code{stat_ellipse(group = 1)}.
#' The original \code{ellipses} parameter (if \code{TRUE}) draws separate
#' ellipses per group defined by the aesthetic mapping; both types of ellipses
#' can appear simultaneously.
#'
#' @return Invisibly, a \code{ggplot} object. Its \code{data} element contains
#'   the first two PCs, the original \code{pheno.df} (if supplied), and if
#'   \code{outlier.pval} is used, additional columns \code{mahal_dist2} and
#'   \code{outlier} (logical).
#' @export
ezpca_outlier <- function(object, pheno.df = NULL, name = "pca", alpha = 1,
                          all.size = NULL, facet = NULL, title = NULL,
                          subtitle = NULL, rm.leg.title = FALSE,
                          labels = FALSE, manual.color = NULL,
                          manual.shape = NULL, ellipses = FALSE,
                          outlier.pval = NULL, label.outliers = TRUE,
                          plot = TRUE, ...) {
  stopifnot(limma::isNumeric(object),
            nrow(object[rowSums(is.na(object)) == 0, ]) > 0,
            ncol(object) > 0,
            !is.null(colnames(object)),
            is.logical(plot))

  # If outlier detection is requested, disable global point labeling
  if (!is.null(outlier.pval)) {
    if (isTRUE(labels)) {
      message("'labels = TRUE' is overridden because outlier detection is active.")
      labels <- FALSE
    }
    stopifnot(is.numeric(outlier.pval), length(outlier.pval) == 1,
              outlier.pval > 0, outlier.pval < 1)
  }

  # PCA on samples (columns), removing rows with NA
  pca <- stats::prcomp(t(object[rowSums(is.na(object)) == 0, ]), scale. = FALSE)
  pve <- signif(summary(pca)$importance["Proportion of Variance", 1:2] * 100, 2)

  # Build data frame with first two PCs and phenotype info
  if (!is.null(pheno.df)) {
    stopifnot(ncol(object) == nrow(pheno.df),
              colnames(object) == rownames(pheno.df))
    dat <- data.frame(pca$x[rownames(pheno.df), 1:2],
                      pheno.df,
                      check.names = FALSE)
    dat$row_names <- rownames(pheno.df)
  } else {
    dat <- data.frame(pca$x[colnames(object), 1:2],
                      check.names = FALSE)
    dat$row_names <- colnames(object)
  }

  # Outlier detection based on Mahalanobis distance (global)
  if (!is.null(outlier.pval)) {
    var1 <- pca$sdev[1]^2
    var2 <- pca$sdev[2]^2
    dat$mahal_dist2 <- (dat$PC1^2) / var1 + (dat$PC2^2) / var2
    threshold <- stats::qchisq(1 - outlier.pval, df = 2)
    dat$outlier <- dat$mahal_dist2 > threshold
  }

  # Adjust PDF width based on longest label in mapped aesthetics (as in original)
  dots <- list(...)
  if (is.null(names(dots))) {
    n <- 0
  } else {
    chars <- vector("list", 2 * length(dots))
    for (i in seq_along(dots)) {
      chars[[2 * i]] <- dots[[i]]
      chars[[2 * i - 1]] <- ifelse(dots[[i]] %in% colnames(dat),
                                   yes = as.character(dat[, dots[[i]]]),
                                   no = "")
    }
    n <- max(nchar(unlist(chars)), na.rm = TRUE)
  }
  width <- 7 + n / 12
  if (!is.na(name)) {
    grDevices::pdf(paste0(name, ".pdf"), width = width, height = 7)
    on.exit(grDevices::dev.off())
  }

  # Base plot
  qp <- ggplot2::ggplot(dat, mapping = ggplot2::aes_string(x = "PC1", y = "PC2", ...)) +
    ggplot2::theme_bw()
  if (!is.null(all.size)) {
    qp <- qp + ggplot2::geom_point(size = all.size, alpha = alpha)
  } else {
    qp <- qp + ggplot2::geom_point(alpha = alpha)
  }
  if (!is.null(facet)) {
    qp <- qp + ggplot2::facet_grid(facet)
  }
  qp <- qp + ggplot2::xlab(paste0("PC1 (", pve[1], "%)")) +
    ggplot2::ylab(paste0("PC2 (", pve[2], "%)"))
  if (rm.leg.title) {
    qp <- qp + ggplot2::theme(legend.title = ggplot2::element_blank())
  }
  if (!is.null(title)) {
    qp <- qp + ggplot2::ggtitle(label = title, subtitle = subtitle)
  }

  # Original ellipses per group (if requested)
  if (ellipses) {
    qp <- qp + ggplot2::stat_ellipse()
  }

  # Outlier ellipse (single group, if requested)
  if (!is.null(outlier.pval)) {
    qp <- qp + ggplot2::stat_ellipse(
      mapping = ggplot2::aes(group = 1),
      level = 1 - outlier.pval,
      type = "norm"   # uses chi‑square distribution
    )
  }

  # Manual scales if provided
  if (!is.null(manual.color)) {
    qp <- qp + ggplot2::scale_colour_manual(values = manual.color)
  }
  if (!is.null(manual.shape)) {
    qp <- qp + ggplot2::scale_shape_manual(values = manual.shape)
  }

  # Global point labels (original feature)
  if (labels) {
    qp <- qp + ggrepel::geom_text_repel(
      data = dat,
      mapping = ggplot2::aes_string(x = "PC1", y = "PC2", label = "row_names"),
      size = 2, vjust = -0.7, show.legend = FALSE
    )
  }

  # Outlier labels (if requested and outliers exist)
  if (!is.null(outlier.pval) && isTRUE(label.outliers) && any(dat$outlier)) {
    qp <- qp + ggrepel::geom_text_repel(
      data = subset(dat, outlier),
      mapping = ggplot2::aes_string(x = "PC1", y = "PC2", label = "row_names"),
      size = 2, vjust = -0.7, show.legend = FALSE
    )
  }

  if (plot) {
    graphics::plot(qp)
  }
  invisible(qp)
}
