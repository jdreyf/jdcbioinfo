#' Find markers by DESeq2's pairwise wald tests
#'
#' Perform DESeq2's pairwise wald tests, find specifically up or down-regulated markers for each group.
#'
#' @param direction Either "up" or "down" or both
#' @inheritParams deseq2_contrasts
#' @inheritParams ezlimma::limma_contrasts
#' @return Data frame.
#' @export

deseq2_find_all_markers <- function(dds, grp, direction=c("up", "down"), design=NULL, add.means=!is.null(grp),
                                adjust.method="BH", ncore=1, shrunken=TRUE){

  stopifnot(ncol(dds)==length(grp), colnames(dds)==names(grp), length(unique(grp))>1)

  # make contrast
  groups <- unique(sort(grp))
  comb <- utils::combn(groups, 2)
  contrasts.v <- character(0)
  for(i in 1:ncol(comb)){
    contrasts.v[paste0(comb[2, i], "_vs_", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  if (!is.null(design)) {
    design(dds) <- design
  }
  mtt <- deseq2_contrasts(dds, grp=grp, contrast.v=contrasts.v, add.means=FALSE, cols=c("stat","pvalue", "padj", "log2FoldChange"),
                          ncore=ncore, shrunken=shrunken)
  mtt <- mtt[, grep("\\.stat$", colnames(mtt))]
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev, ".stat")
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  resAll <- list()
  for (d in direction) {
    res <- list()
    score_fn <- switch(d, up=min, down=max)

    for(i in seq_along(groups)){
      nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
      mtt_tmp <- mtt[, paste0(nms, ".stat"), drop=FALSE]
      score <- apply(mtt_tmp, 1, score_fn, na.rm=TRUE)
      score[is.infinite(score)] <- NA # for min/max of all NAs

      if(d=="up"){
        pval <- apply(mtt_tmp, MARGIN = 1, FUN = function(scores) prod(1 - stats::pnorm(scores)))
      }else if(d=="down"){
        pval <-  apply(mtt_tmp, MARGIN = 1, FUN = function(scores) prod(stats::pnorm(scores)))
      }

      fdr <- stats::p.adjust(pval, method=adjust.method)
      res_tmp <- data.frame(score=score, p=pval, FDR=fdr)
      colnames(res_tmp) <- paste(groups[i], d, colnames(res_tmp), sep=".")
      res[[i]] <- res_tmp
    }
    res <- Reduce(cbind, res)
    res <- res[rownames(dds), ]
    resAll[[d]] <- res
  }

  resAll <- Reduce(cbind, resAll)
  if(add.means){
    mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
    mat_avg <- sapply(groups, function(g) rowMeans(mat[, grp==g, drop=FALSE]))
    colnames(mat_avg) <- paste0(groups, ".avg")
    resAll <- cbind(mat_avg[rownames(resAll), ], resAll)
  }
  return(resAll)
}
