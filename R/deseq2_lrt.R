#' Apply DESeq2's DESeq, results, and lfcShrink to one or more contrasts, and return a data.frame
#'
#' Apply \pkg{DESeq2}'s \code{DESeq}, \code{results},and \code{lfcShrink} to one or more contrasts, and return a data.frame
#' @param dds DESeqDataSet object.
#' @param reduced model matrix of the reduced model
#' @param prefix the prefix for the test result columns
#' @param ncore Number of cores to use.
#' @inheritParams ezlimma::limma_contrasts
#' @details  \code{grp} isn't needed if \code{add.means} is \code{FALSE}.
#' @export

deseq2_lrt <- function(dds, grp=NULL, reduced=design(dds)[,1], prefix = "Group", add.means=!is.null(grp),
                       cols=c("stat", "pvalue", "padj"), ncore=1){

    if (add.means) { stopifnot(ncol(dds)==length(grp), colnames(dds)==names(grp)) }

    if (.Platform$OS.type=="windows") {
      bp <- BiocParallel::SnowParam(workers=ncore, type="SOCK")
    } else {
      bp <- BiocParallel::MulticoreParam(workers=ncore)
    }

    BiocParallel::register(BiocParallel::bpstart(bp))

    # Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
    dds <- DESeq2::DESeq(dds, test="LRT", reduced=reduced, parallel=TRUE, BPPARAM=bp)

    # results
    tt <- DESeq2::results(dds, cooksCutoff=FALSE, independentFiltering=FALSE, test="LRT", parallel=TRUE, BPPARAM=bp)

    tt <- as.data.frame(tt)[, cols]
    colnames(tt) <- gsub("pvalue", "p", gsub("padj", "FDR", colnames(tt)))
    colnames(tt) <- paste

    BiocParallel::bpstop(bp)
    tt <- tt[order(tt[, paste0(prefix, ".p")]), ](prefix, colnames(tt), sep=".")

    if (add.means) {
      mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
      groups <- unique(sort(grp))
      grp.means <- vapply(groups, function(g) rowMeans(mat[, grp==g, drop = FALSE]), numeric(nrow(mat)))
      colnames(grp.means) <- paste(groups, "avg", sep=".")
      tt <- cbind(grp.means[rownames(tt),], tt)
    }
    return(tt)
}
