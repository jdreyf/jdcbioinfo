#' Plot logit missing probability by average log abundance
#'
#' Plot logit missing probability by average log abundance.
#'
#' @param mat Numeric matrix of log abundance.
#' @param name file names of the plot.
#' @param main title of the plot.
#' @param alpha Transparency of color.
#' @param offset offset of the number of missing values when calculating the proportion of NAs with
#'  `(rowSums(is.na(mat))+offset)/(ncol(mat)+offset)`, where need to avoid the proportion of NAs being zero
#'  so that we can calculate its logit.
#' @return NULL.
#' @export

plot_missing <- function(mat, name=NA, main="", alpha=1, offset=0.5){
  stopifnot(!is.null(rownames(mat)))

  na.prop <- (rowSums(is.na(mat))+offset)/(ncol(mat)+offset)
  na.prop <- na.prop[na.prop != 1]
  na.logit <- log2(na.prop/(1 - na.prop))
  avg.log <- rowMeans(mat[names(na.prop), ], na.rm=TRUE)

  lm1 <- stats::lm(na.logit ~ avg.log)
  slope <- signif(stats::coef(summary(lm1))["avg.log", "Estimate"], 3)
  pval <- signif(stats::coef(summary(lm1))["avg.log", "Pr(>|t|)"], 3)

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"))
    on.exit(grDevices::dev.off())
  }

  ord <- ifelse(slope > 0, "topleft", "topright")
  graphics::plot(x=avg.log, y=jitter(na.logit, 2), xlab="Average log abundance", ylab="Logit missingness probability",
                 pch=16, col=scales::alpha("black", alpha=alpha), main=main)
  graphics::abline(lm1, col="red")
  graphics::legend(ord, legend=paste0("slope = ", slope, "\npval = ", pval), bty="n")
}
