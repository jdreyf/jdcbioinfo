#' Find markers for each group by limma's pairwise moderated t-tests
#'
#' Perform limma's pairwise moderated t-tests, find specifically up or down-regulated markers for each group.
#'
#' @param direction Either "up" or "down" or both
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams rankprod
#' @return Data frame.
#' @export

limma_find_all_markers <- function(object, grp, direction= c("up", "down"), design=NULL, add.means=!is.null(grp),
                                   adjust.method="BH", weights=NA, trend=FALSE, block=NULL, correlation=NULL,
                                   treat.lfc=NULL, moderated=TRUE){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp), length(unique(grp))>1, nsim<=1e7-2)

  # make contrast
  groups <- unique(sort(grp))
  comb <- utils::combn(groups, 2)
  contrasts.v <- character(0)
  for(i in 1:ncol(comb)){
    contrasts.v[paste0(comb[2, i], "_vs_", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  mtt <- ezlimma::limma_contrasts(object, grp=grp, contrast.v=contrasts.v, design=design, weights=weights, trend=trend, block=block,
                                  correlation=correlation, adjust.method=adjust.method, add.means=FALSE, treat.lfc=treat.lfc, moderated=moderated,
                                  check.names=TRUE, cols=c("P.Value", "logCF"))
  mtt <- multi_pval2z(mtt)
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev)
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  resAll <- list()
  for (d in direction) {
    res <- list()
    score_fn <- switch(d, up=min, down=max)

    for(i in seq_along(groups)){
      nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
      mtt_tmp <- mtt[, paste0(nms)]
      score <- apply(mtt_tmp, 1, score_fn, na.rm=TRUE)
      score[is.infinite(score)] <- NA # for min/max of all NAs

      n <- length(groups)-1
      if(d=="up"){
        pval <- (1 - stats::pnorm(score))^n
      }else if(d=="down"){
        pval <- stats::pnorm(score)^n
      }

      fdr <- stats::p.adjust(pval, method=adjust.method)
      res_tmp <- data.frame(score=score, p=pval, FDR=fdr)
      colnames(res_tmp) <- paste(groups[i], direction, colnames(res_tmp), sep=".")
      res[[i]] <- res_tmp
    }
    res <- Reduce(cbind, res)
    res <- res[rownames(object), ]
    resAll[[d]] <- res
  }

  resAll <- Reduce(cbind, resAll)
  if(add.means){
    mat_avg <- sapply(groups, function(g) rowMeans(object[, grp==g]))
    colnames(mat_avg) <- paste0(groups, ".avg")
    resAll <- cbind(mat_avg[rownames(resAll), ], resAll)
  }
  return(resAll)
}
