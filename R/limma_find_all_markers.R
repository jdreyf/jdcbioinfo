#' Find markers for each group by limma's pairwise moderated t-tests
#'
#' Perform limma's pairwise moderated t-tests, find specifically up and down-regulated markers for each group.
#'
#' @inheritParams ezlimma::limma_contrasts
#' @return Data frame.
#' @export

limma_find_all_markers <- function(object, grp, design=NULL, add.means=!is.null(grp),
                                   adjust.method="BH", weights=NA, trend=FALSE, block=NULL, correlation=NULL,
                                   treat.lfc=NULL, moderated=TRUE){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp), length(unique(grp))>1)

  # make contrast
  groups <- unique(sort(grp))
  comb <- utils::combn(groups, 2)
  contrasts.v <- character(0)
  for(i in 1:ncol(comb)){
    contrasts.v[paste0(comb[2, i], "_vs_", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  mtt <- ezlimma::limma_contrasts(object, grp=grp, contrast.v=contrasts.v, design=design, weights=weights, trend=trend, block=block,
                                  correlation=correlation, adjust.method=adjust.method, add.means=FALSE, treat.lfc=treat.lfc, moderated=moderated,
                                  check.names=TRUE, cols=c("P.Value", "logFC"))
  mtt <- multi_pval2z(mtt)
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev)
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  n <- length(groups)-1
  res <- list()
  for(i in seq_along(groups)){
    nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
    mtt_tmp <- mtt[, paste0(nms)]

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
  res <- res[rownames(object), ]

  if(add.means){
    mat_avg <- sapply(groups, function(g) rowMeans(object[, grp==g, drop = FALSE]))
    colnames(mat_avg) <- paste0(groups, ".avg")
    res <- cbind(mat_avg[rownames(res), ], res)
  }
  return(res)
}
