#' Find markers for each group by limma's pairwise moderated t-tests over a logFC cutoff
#'
#' Perform limma's pairwise moderated t-tests, find specifically up or down-regulated markers for each group over a logFC cutoff.
#'
#' @param direction Either "up" or "down" or both
#' @inheritParams ezlimma::limma_contrasts
#' @inheritParams rankprod
#' @return Data frame.
#' @export

limma_treat_find_all_markers <- function(object, grp, direction= c("up", "down"), treat.lfc=log2(1.2), design=NULL,
                                         add.means=!is.null(grp),
                                         adjust.method="BH", weights=NA, trend=FALSE, block=NULL, correlation=NULL,
                                         moderated=TRUE){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp), length(unique(grp))>1)

  # make model
  if (is.null(design)){
    design <- stats::model.matrix(~0+grp)
    colnames(design) <- sub('grp', '', colnames(design), fixed=TRUE)
  }

  # lmFit
  if (!is.na(weights)){
    if (!is.matrix(object) && !is.null(object$weights)){ cat('object$weights are being ignored\n') }
    fit <- limma::lmFit(object, design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- limma::lmFit(object, design, block = block, correlation = correlation)
  }

  # make contrast
  groups <- unique(sort(grp))
  comb <- utils::combn(groups, 2)
  contrasts.v <- character(0)
  for(i in 1:ncol(comb)){
    contrasts.v[paste0(comb[2, i], "_vs_", comb[1, i])] <- paste0(comb[2, i], " - ", comb[1, i])
  }

  contr.mat <- limma::makeContrasts(contrasts=contrasts.v, levels=design)

  # contrasts.fit & treat
  fit2 <- limma::contrasts.fit(fit, contr.mat)
  fit2 <- limma::treat(fit2, lfc=treat.lfc, trend=trend)

  # topTable
  for (i in seq_along(contrasts.v)) {
    tt <- limma::topTreat(fit2, number=Inf, coef=contrasts.v[i])
    tt <- tt[, c("P.Value", "logFC")]
    colnames(tt) <- paste(names(contrasts.v)[i], c("p", "logFC"), sep = ".")

    if(i == 1) {
      mtt <- tt
    } else {
      mtt <- cbind(mtt, tt[rownames(mtt), ])
    }
  }

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
      mtt_tmp <- mtt[, nms, drop=FALSE]
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
