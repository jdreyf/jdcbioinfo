#' #' Apply edgeR's glmQLFit or glmFit, glmQLFTest, glmLRT, or glmWeightedF to one or more contrasts, perform F-test for all contrasts,  and return a list
#'
#' #'Apply \pkg{edgeR}'s \code{glmQLFit}, or \code{glmFit},  \code{glmQLFTest}, \code{glmLRT}, or \code{glmWeightedF} to one or more contrasts,
#' perform F-test for all contrasts, and return a list
#'
#' @inheritParams edger_contrasts
#' @inheritParams limma_ftest_contrasts
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means}
#' is \code{FALSE}.
#' @export

edger_ftest_contrasts <- function(dge, grp=NULL, contrast.v, add.means=!is.null(grp), design=NULL, prefix="", cols=c("F", "t", "PValue", "FDR"),
                                  test="QLFT", plot=FALSE){

    if (is.null(design)|add.means) { stopifnot(ncol(dge)==length(grp), colnames(dge)==names(grp)) }
    stopifnot(test %in% c("QLFT", "LRT",  "WeightedFT"))

    if (is.null(design)){
        design <- stats::model.matrix(~0+grp)
        colnames(design) <- sub("grp", "", colnames(design), fixed=TRUE)
    }

    dge <- edgeR::estimateDisp(dge, design=design)
    if (plot) edgeR::plotBCV(dge)

    if (test=="QLFT") {
        fit <- edgeR::glmQLFit(dge, design=design)
    } else if (test=="LRT") {
        fit <- edgeR::glmFit(dge, design=design)
    } else if (test=="WeightedFT") {
        fit <- edgeR::glmFit(dge, design=design)
        base_mean <- unname(rowMeans(sweep(dge$counts, 2, dge$samples$norm.factors, FUN="*")))
    }
    if (plot && test=="QLFT") edgeR::plotQLDisp(fit)

    contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=design)

    if (test=="QLFT") {
        tres <- edgeR::glmQLFTest(fit, contrast=contr.mat)
    } else if (test=="LRT") {
        tres <- edgeR::glmLRT(fit, contrast=contr.mat)
    } else if (test=="WeightedFT") {
        tres <- zinbwave::glmWeightedF(fit, contrast=contr.mat, filter=base_mean)
    }
    if (plot) limma::plotMD(tres)

    mtt <- edgeR::topTags(tres, n=nrow(dge), adjust.method="BH")
    patt <- paste(paste0('^', cols, '$'), collapse='|')
    mtt <- as.data.frame(mtt[, grep(patt, colnames(mtt))])
    colnames(mtt) <-gsub("PValue", "p", colnames(mtt))
    if (prefix!=""){ colnames(mtt) <- paste(prefix, colnames(mtt), sep=".") }

    logcpm <- edgeR::cpm(dge, normalized.lib.sizes=TRUE, log=TRUE, prior.count=fit$prior.count)

    if (add.means) {
        groups <- unique(sort(grp))
        grp.means <- vapply(groups, function(g) rowMeans(logcpm[, grp==g]), numeric(nrow(logcpm)))
        colnames(grp.means) <- paste(groups, "avg.logcpm", sep=".")
        mtt <- cbind(grp.means[rownames(mtt),], mtt)
    }

    return(list(mtt=mtt, logcpm=logcpm, dge_obj=dge))
}