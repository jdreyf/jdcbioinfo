#' Find markers by DESeq2's pairwise wald tests
#'
#' Perform DESeq2's pairwise wald tests, find specifically up or down-regulated markers for each group.
#'
#' @param direction Either "up" or "down".
#' @inheritParams deseq2_contrasts
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams rankprod
#' @return Data frame.
#' @export

deseq2_find_markers <- function(dds, grp, direction= c("up", "down"), nsim=1e7-2, seed=100, design=NULL,
                                add.means=!is.null(grp),
                                adjust.method="BH", ncore=1, shrunken=TRUE){

  warning("deseq2_find_markers is depreciated, please use deseq2_find_all_markers")
  stopifnot(ncol(dds)==length(grp), colnames(dds)==names(grp), length(unique(grp))>1, nsim<=1e7-2)

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
    set.seed(seed)
    for(i in seq_along(groups)){
      nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
      mtt_tmp <- mtt[, paste0(nms, ".stat")]
      score_stat <- apply(mtt_tmp, 1, score_fn, na.rm=TRUE)
      score_stat[is.infinite(score_stat)] <- NA # for min/max of all NAs

      mtt_sim <- apply(mtt_tmp, 2, function(v) sample(v, nsim, replace=TRUE))
      score_stat_sim <- apply(mtt_sim, 1, score_fn, na.rm=TRUE)
      score_stat_sim[is.infinite(score_stat_sim)] <- NA
      Fn <- stats::ecdf(c(score_stat_sim, Inf, -Inf))

      if(d=="up"){
        pval <- 1 - Fn(score_stat)
      }else if(d=="down"){
        pval <- Fn(score_stat)
      }

      fdr <- stats::p.adjust(pval, method=adjust.method)
      res_tmp <- data.frame(score=score_stat, p=pval, FDR=fdr)
      colnames(res_tmp) <- paste(groups[i], direction, colnames(res_tmp), sep=".")
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
