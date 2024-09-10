#' Apply DESeq2's DESeq, results, and lfcShrink to one or more contrasts, and return a data.frame
#'
#' Apply \pkg{DESeq2}'s \code{DESeq}, \code{results},and \code{lfcShrink} to one or more contrasts, and return a data.frame
#' @param dds DESeqDataSet object.
#' @param ncore Number of cores to use.
#' @param shrunken logical, whether to shrink the log fold-change
#' @inheritParams ezlimma::limma_contrasts
#' @details  \code{grp} isn't needed if \code{add.means} is \code{FALSE}.
#' @export

deseq2_contrasts <- function(dds, grp=NULL, contrast.v, add.means=!is.null(grp), cols=c("pvalue", "padj", "log2FoldChange"), ncore=1, shrunken=TRUE){

    if (add.means) { stopifnot(ncol(dds)==length(grp), colnames(dds)==names(grp)) }

    if (.Platform$OS.type=="windows") {
      bp <- BiocParallel::SnowParam(workers=ncore, type="SOCK")
    } else {
      bp <- BiocParallel::MulticoreParam(workers=ncore)
    }

    BiocParallel::register(BiocParallel::bpstart(bp))

    # Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
    dds <- DESeq2::DESeq(dds, test="Wald", parallel=TRUE, BPPARAM=bp)

    # contrasts
    contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=BiocGenerics::design(dds))
    colnames(contr.mat) <- names(contrast.v)

    for (i in 1:ncol(contr.mat)) {
      tt <- DESeq2::results(dds, contrast=contr.mat[, i], cooksCutoff=FALSE, independentFiltering=FALSE,
                            test="Wald", parallel=TRUE, BPPARAM=bp)
      if (shrunken) {
        shrunkenTT <- DESeq2::lfcShrink(dds, contrast=contr.mat[, i], res=tt, type="ashr", parallel=TRUE, BPPARAM=bp)
      } else {
        shrunkenTT <- tt
      }

      shrunkenTT <- as.data.frame(shrunkenTT)[, cols]
      colnames(shrunkenTT) <- gsub("pvalue", "p", gsub("padj", "FDR", gsub("log2FoldChange", "logFC", colnames(shrunkenTT))))
      shrunkenTT$FC <- sign(shrunkenTT$logFC)*2^abs(shrunkenTT$logFC)
      colnames(shrunkenTT) <- paste(colnames(contr.mat)[i], colnames(shrunkenTT), sep=".")

      if (i == 1) {
        mtt <- shrunkenTT
      } else {
        mtt <- cbind(mtt, shrunkenTT[rownames(mtt), ])
      }
    }
    BiocParallel::bpstop(bp)
    mtt <- mtt[order(ezlimma::combine_pvalues(mtt)), ]

    mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))

    if (add.means) {
        groups <- unique(sort(grp))
        grp.means <- vapply(groups, function(g) rowMeans(mat[, grp==g, drop = FALSE]), numeric(nrow(mat)))
        colnames(grp.means) <- paste(groups, "avg", sep=".")
        mtt <- cbind(grp.means[rownames(mtt),], mtt)
    }
    return(mtt)
}
