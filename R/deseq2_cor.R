#' Apply DESeq2's DESeq, results, and lfcShrink to one or more contrasts, and return a data.frame
#'
#' Apply \pkg{DESeq2}'s \code{DESeq}, \code{results},and \code{lfcShrink} to one or more contrasts, and return a data.frame
#' @param dds DESeqDataSet object.
#' @param name Name of the individual effect (coefficient) to test
#' @param ncore Number of cores to use.
#' @param shrunken logical, whether to shrink the slope
#' @inheritParams ezlimma::limma_contrasts
#' @export

deseq2_cor <- function(dds, name = colnames(design(dds))[2], cols=c("pvalue", "padj", "log2FoldChange"), ncore=1, shrunken=TRUE){

    if (.Platform$OS.type=="windows") {
      bp <- BiocParallel::SnowParam(workers=ncore, type="SOCK")
    } else {
      bp <- BiocParallel::MulticoreParam(workers=ncore)
    }

    BiocParallel::register(BiocParallel::bpstart(bp))

    # Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
    dds <- DESeq2::DESeq(dds, test="Wald", parallel=TRUE, BPPARAM=bp)

    tt <- DESeq2::results(dds, name=name, cooksCutoff=FALSE, independentFiltering=FALSE,
                          test="Wald", parallel=TRUE, BPPARAM=bp)
    if (shrunken) {
      shrunkenTT <- DESeq2::lfcShrink(dds, coef=name, res=tt, type="ashr", parallel=TRUE, BPPARAM=bp)
    } else {
      shrunkenTT <- tt
    }

    shrunkenTT <- as.data.frame(shrunkenTT)[, cols]
    colnames(shrunkenTT) <- gsub("pvalue", "p", gsub("padj", "FDR", gsub("log2FoldChange", "slope", colnames(shrunkenTT))))
    colnames(shrunkenTT) <- paste(name, colnames(shrunkenTT), sep=".")

    BiocParallel::bpstop(bp)

    return(shrunkenTT)
}
