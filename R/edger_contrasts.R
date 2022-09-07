#' Apply edgeR's glmQLFit or glmFit, glmQLFTest, glmLRT, or glmWeightedF to one or more contrasts, and return a list
#'
#' Apply \pkg{edgeR}'s \code{glmQLFit}, or \code{glmFit},  \code{glmQLFTest}, \code{glmLRT}, or \code{glmWeightedF} to one or more contrasts, and return
#' a list
#' @param dge DGEList object.
#' @param test QLFT, LRT or WeightedFT.
#' @param plot Logical indicate if to plot the QC metrics.
#' @inheritParams ezlimma::limma_contrasts
#' @details If \code{design} is \code{NULL} and \code{grp} is given, design will be calculated as
#' \code{model.matrix(~0+grp)}. However, \code{grp} isn't needed if \code{design} is provided & \code{add.means}
#' is \code{FALSE}.
#' @export

edger_contrasts <- function(dge, grp=NULL, contrast.v, add.means=!is.null(grp), design=NULL, cols=c("PValue", "FDR", "logFC"),
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
    for (i in seq_along(contrast.v)) {
        if (test=="QLFT") {
            tres <- edgeR::glmQLFTest(fit, contrast=contr.mat[, i])
        } else if (test=="LRT") {
            tres <- edgeR::glmLRT(fit, contrast=contr.mat[, i])
        } else if (test=="WeightedFT") {
            if (!requireNamespace("zinbwave", quietly = TRUE)) stop("Package \"zinbwave\" must be installed to use this function.", call. = FALSE)
            tres <- zinbwave::glmWeightedF(fit, contrast=contr.mat[, i], filter=base_mean)
        }
        if (plot) limma::plotMD(tres)

        tt <- edgeR::topTags(tres, n=nrow(dge), adjust.method="BH")
        tt <- as.data.frame(tt[, cols])
        tt$FC <- sign(tt$logFC) * 2^abs(tt$logFC)
        colnames(tt) <- paste(names(contrast.v)[i], gsub("PValue", "p", colnames(tt)), sep=".")
        if (i==1) { mtt <- tt }
        else { mtt <- cbind(mtt, tt[rownames(mtt), ]) }
    }

    mtt <- mtt[order(ezlimma::combine_pvalues(mtt)), ]
    logcpm <- edgeR::cpm(dge, normalized.lib.sizes=TRUE, log=TRUE, prior.count=fit$prior.count)

    if (add.means) {
        groups <- unique(sort(grp))
        grp.means <- vapply(groups, function(g) rowMeans(logcpm[, grp==g]), numeric(nrow(logcpm)))
        colnames(grp.means) <- paste(groups, "avg", sep=".")
        mtt <- cbind(grp.means[rownames(mtt),], mtt)
    }

    return(list(mtt=mtt, logcpm=logcpm, dge_obj=dge))
}
