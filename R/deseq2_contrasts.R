#' Apply DESeq2's DESeq, results, and lfcShrink to one or more contrasts, and return a data.frame
#'
#' Apply \pkg{DESeq2}'s \code{DESeq}, \code{results},and \code{lfcShrink} to one or more contrasts, and return a data.frame
#' @param dds DESeqDataSet object.
#' @param ncore Number of cores to use.
#' @param shrunken logical, whether to shrink the log fold-change.
#' @param add_cooks logical, whether to include Cook's cutoff. \code{DESeq2} default is \code{TRUE}.
#' @param lfc_shrink_type character string, method for \code{lfcShrink}. It is recommended to not use method \code{"normal"}.
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams DESeq2::results
#' @details \code{grp} isn't needed if \code{add.means} is \code{FALSE}. To match \code{DESeq2} default, set \code{add_cooks = independentFiltering = TRUE}.
#' This function's default values are chosen to be back-compatible with previous versions.
#' @export

deseq2_contrasts <- function(dds, grp=NULL, contrast.v, add.means=!is.null(grp),
                             cols=c("pvalue", "padj", "log2FoldChange"),
                             shrunken=FALSE,
                             add_cooks=FALSE,
                             independentFiltering=FALSE,
                             parallel=TRUE,
                             ncore=1,
                             alpha=0.1,
                             lfc_shrink_type=c("ashr", "apeglm", "normal")){

    lfc_shrink_type <- match.arg(lfc_shrink_type)
    if (add.means) { stopifnot(ncol(dds)==length(grp), colnames(dds)==names(grp)) }

    if (parallel){
      if (.Platform$OS.type=="windows") {
        bp <- BiocParallel::SnowParam(workers=ncore, type="SOCK")
      } else {
        bp <- BiocParallel::MulticoreParam(workers=ncore)
      }
      BiocParallel::register(BiocParallel::bpstart(bp))
    }

    # Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
    if (parallel){
      dds <- DESeq2::DESeq(dds, test="Wald", parallel=TRUE, BPPARAM=bp)
    } else {
      dds <- DESeq2::DESeq(dds, test="Wald")
    }

    # contrasts
    contr.mat <- limma::makeContrasts(contrasts=contrast.v, levels=BiocGenerics::design(dds))
    colnames(contr.mat) <- names(contrast.v)

    for (i in 1:ncol(contr.mat)) {
      args_lst <- list(object = dds, contrast=contr.mat[, i], independentFiltering=independentFiltering,
                       test="Wald", alpha=alpha)
      if (parallel) args_lst <- c(args_lst, parallel=TRUE, BPPARAM=bp)
      if (!add_cooks) args_lst <- c(args_lst, cooksCutoff=FALSE)

      tt <- do.call(DESeq2::results, args = args_lst)

      if (shrunken & parallel) {
        shrunkenTT <- DESeq2::lfcShrink(dds, contrast=contr.mat[, i], res=tt, type=lfc_shrink_type, parallel=TRUE, BPPARAM=bp)
      } else if (shrunken & !parallel){
        shrunkenTT <- DESeq2::lfcShrink(dds, contrast=contr.mat[, i], res=tt, type=lfc_shrink_type)
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
    if (parallel) BiocParallel::bpstop(bp)
    mtt <- mtt[order(ezlimma::combine_pvalues(mtt)), ]

    if (add.means) {
      mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
      groups <- unique(sort(grp))
      grp.means <- vapply(groups, function(g) rowMeans(mat[, grp==g, drop = FALSE]), numeric(nrow(mat)))
      colnames(grp.means) <- paste(groups, "avg", sep=".")
      mtt <- cbind(grp.means[rownames(mtt),], mtt)
    }
    return(mtt)
}
