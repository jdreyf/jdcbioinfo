#' Find markers by DESeq2's pairwise wald tests
#'
#' Perform DESeq2's pairwise wald tests, find specifically up and down-regulated markers for each group.
#'
#' @inheritParams deseq2_contrasts
#' @inheritParams ezlimma::limma_contrasts
#' @return Data frame.
#' @export

deseq2_find_all_markers <- function(dds, grp, design=NULL, add.means=!is.null(grp),
                                adjust.method="BH", ncore=1){

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
                          ncore=ncore, shrunken=FALSE)
  mtt <- mtt[, grep("\\.stat$", colnames(mtt))]
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev, ".stat")
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  n <- length(groups)-1
  res <- list()
  for(i in seq_along(groups)){
    nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
    mtt_tmp <- mtt[, paste0(nms, ".stat")]

    up_score <- apply(mtt_tmp, 1, min, na.rm=TRUE)
    up_score[is.infinite(up_score)] <- NA # for min/max of all NAs
    dn_score <- apply(mtt_tmp, 1, max, na.rm=TRUE)
    dn_score[is.infinite(dn_score)] <- NA

    score_tab <- cbind(up_score, dn_score)
    score_rev_tab <- cbind(up_score, -1*dn_score)
    pidx <- apply(score_rev_tab, MARGIN = 1, FUN = which.max)

    score <- unlist(Map(function(rnum, idx) {
      score_tab[rnum, idx]
    }, 1:nrow(score_tab), pidx))
    score[score == numeric(0)] <- NA

    score_rev <- unlist(Map(function(rnum, idx) {
      score_rev_tab[rnum, idx]
    }, 1:nrow(score_rev_tab), pidx))
    score_rev[score_rev == numeric(0)] <- NA

    pval <- (stats::pnorm(-1*score_rev))^n
    fdr <- stats::p.adjust(pval, method=adjust.method)
    res_tmp <- data.frame(score=score, p=pval, FDR=fdr)
    colnames(res_tmp) <- paste(groups[i], colnames(res_tmp), sep=".")
    res[[i]] <- res_tmp
  }
  res <- Reduce(cbind, res)
  rownames(res) <- rownames(mtt)
  res <- res[rownames(object), ]

  if(add.means){
    mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
    mat_avg <- sapply(groups, function(g) rowMeans(mat[, grp==g, drop=FALSE]))
    colnames(mat_avg) <- paste0(groups, ".avg")
    res <- cbind(mat_avg[rownames(res), ], res)
  }
  return(res)
}
