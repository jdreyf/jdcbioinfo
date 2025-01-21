#' Find markers by limma's pairwise moderated t-tests
#'
#' Perform limma's pairwise moderated t-tests, find sepcifically up or down-regulated markers for each group.
#'
#' @param direction Either "up" or "down".
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams rankprod
#' @return Data frame.
#' @export

limma_find_markers <- function(object, grp, direction= c("up", "down"), nsim=1e7-2, seed=100, design=NULL,
                               add.means=!is.null(grp),
                               adjust.method="BH", weights=NA, trend=FALSE, block=NULL, correlation=NULL,
                               moderated=TRUE){

  warning("limma_find_markers is depreciated, please use limma_find_all_markers", call. = FALSE)
  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp), length(unique(grp))>1, nsim<=1e7-2)

  # make contrast
  groups <- unique(sort(grp))
  comb <- utils::combn(groups, 2)
  contrasts.v <- character(0)
  for(i in 1:ncol(comb)){
    contrasts.v[paste0(comb[2, i], "_vs_", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  mtt <- ezlimma::limma_contrasts(object, grp=grp, contrast.v=contrasts.v, design=design, weights=weights, trend=trend, block=block,
                                  correlation=correlation, adjust.method=adjust.method, add.means=FALSE, treat.lfc=NULL, moderated=moderated,
                                  check.names=TRUE, cols=c("t", "P.Value"))
  mtt <- mtt[, grep("\\.t$", colnames(mtt))]
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev, '.t')
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  resAll <- list()
  for (d in direction) {
    res <- list()
    score_fn <- switch(d, up=min, down=max)
    set.seed(seed)
    for(i in seq_along(groups)){
      nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
      mtt_tmp <- mtt[, paste0(nms, '.t')]
      score_t <- apply(mtt_tmp, 1, score_fn, na.rm=TRUE)
      score_t[is.infinite(score_t)] <- NA # for min/max of all NAs

      mtt_sim <- apply(mtt_tmp, 2, function(v) sample(v, nsim, replace=TRUE))
      score_t_sim <- apply(mtt_sim, 1, score_fn, na.rm=TRUE)
      score_t_sim[is.infinite(score_t_sim)] <- NA
      Fn <- stats::ecdf(c(score_t_sim, Inf, -Inf))

      if(d=="up"){
        pval <- 1 - Fn(score_t)
      }else if(d=="down"){
        pval <- Fn(score_t)
      }

      fdr <- stats::p.adjust(pval, method=adjust.method)
      res_tmp <- data.frame(score=score_t, p=pval, FDR=fdr)
      colnames(res_tmp) <- paste(groups[i], d, colnames(res_tmp), sep=".")
      res[[i]] <- res_tmp
    }
    res <- Reduce(cbind, res)
    res <- res[rownames(object), ]
    resAll[[d]] <- res
  }

  resAll <- Reduce(cbind, resAll)
  if(add.means){
    mat_avg <-  t(apply(object, 1, FUN=function(v) tapply(v, grp, mean, na.rm=TRUE)))
    colnames(mat_avg) <- paste0(colnames(mat_avg), ".avg")
    resAll <- cbind(mat_avg[rownames(resAll), ], resAll)
  }
  return(resAll)
}
