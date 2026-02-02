#' Find markers by DESeq2's pairwise wald tests
#'
#' Perform DESeq2's pairwise wald tests, find specifically up or down-regulated markers for each group.
#'
#' @param direction Either "up" or "down" or both
#' @param dds A DESeqDataSet object
#' @param grp A named vector of group labels for each column in dds
#' @param design Optional design formula to set for dds
#' @param method Method for p-value calculation, one of "exact", "asymptotic", or "monte_carlo"
#' @param nsim Number of Monte Carlo simulations to perform if method is "monte_carlo"
#' @param seed Random seed for reproducibility
#' @param add.means Whether to add group-wise mean expression values to the output
#' @param adjust.method Method for p-value adjustment, passed to p.adjust
#' @param ncore Number of cores to use
#' @param shrunken Whether to use shrunken log2 fold changes for ranking
#' @return Data frame.
#' @export

deseq2_find_all_markers <- function(dds, grp, direction=c("up", "down"), design=NULL,
                                    method=c("exact", "asymptotic", "monte_carlo"),
                                    nsim=1e7-2, seed=42, add.means=!is.null(grp),
                                    adjust.method="BH", ncore=1, shrunken=FALSE){

  stopifnot(ncol(dds)==length(grp),
            colnames(dds)==names(grp),
            length(unique(grp))>1)
  method <- match.arg(method)
  set.seed(seed)

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
  mtt <- deseq2_contrasts(dds,
                          grp=grp,
                          contrast.v=contrasts.v,
                          add.means=FALSE,
                          cols=c("stat","pvalue", "padj", "log2FoldChange"),
                          ncore=ncore,
                          shrunken=shrunken)
  mtt <- mtt[, grep("\\.stat$", colnames(mtt))]
  mtt_rev <- -1*mtt
  nms_rev <- sapply(1:ncol(comb), function(i) paste0(comb[1, i], "_vs_", comb[2, i]))
  colnames(mtt_rev) <- paste0(nms_rev, ".stat")
  mtt <- cbind(mtt, mtt_rev)
  rm(mtt_rev)

  k <- length(groups) - 1
  if(method=="asymptotic") {
    # Warn about k
    if(k<10){
      warning("k = ", k, ": Asymptotic method is NOT RECOMMENDED for k < 10. Use exact or Monte Carlo instead.")
    }else if(k < 20){
      warning("k = ", k, ": Asymptotic method may be INACCURATE for k < 20. Consider exact method for better accuracy.")
    }else if(k<50){
      warning("k = ", k, ": Asymptotic method is REASONABLE for k >= 20 but may have bias in extreme tails.")
    }

    # Precompute Gumbel parameters
    if(k > 1){
      b_n<-sqrt(2*log(k))-(log(log(k))+log(4*pi))/(2*sqrt(2*log(k)))
      a_n <- 1/sqrt(2*log(k))
    }else{
      b_n <- 0
      a_n <- 1
    }
  }

  res_all <- list()
  for(d in direction){
    res <- list()
    score_fn <- switch(d, up=min, down=max)

    for(i in seq_along(groups)){
      nms <- sapply(setdiff(seq_along(groups), i), function(j) paste0(groups[i], "_vs_", groups[j]))
      mtt_tmp <- mtt[, paste0(nms, ".stat"), drop=FALSE]
      score <- apply(mtt_tmp, 1, score_fn, na.rm=TRUE)
      score[is.infinite(score)] <- NA # for min/max of all NAs

      if(method=="exact"){
        if(d=="up"){
          # For positive minimum: P(min ≥ observed) = (1 - Φ(score))^k
          pval <- (1 - stats::pnorm(score))^k
        }else if(d=="down"){
          # For negative maximum: P(max ≤ observed) = Φ(score)^k
          pval <- (stats::pnorm(score))^k
        }
      }else if(method=="asymptotic"){
        if(k==1){
          # For k=1, it's just standard normal
          if(d=="up") {
            pval <- 1 - stats::pnorm(score)  # Right tail
          }else if(d=="down") {
            pval <- stats::pnorm(score)      # Left tail
          }
        }else{
          if(d=="up") {
            # For positive minimum: use Gumbel for -score
            # P(min ≥ score) = P(max(-X) ≤ -score)
            gumbel_val <- (-score - b_n) / a_n
            pval <- evd::pgumbel(gumbel_val, loc = 0, scale = 1)
          }else if(d=="down") {
            # For negative maximum: use Gumbel directly
            # P(max ≤ score) = F_gumbel(score)
            gumbel_val <- (score - b_n) / a_n
            pval <- evd::pgumbel(gumbel_val, loc = 0, scale = 1)
          }
        }
      }else if(method == "monte_carlo"){
        mtt_sim <- apply(mtt_tmp, 2, function(v) sample(v, nsim, replace=TRUE))
        score_sim <- apply(mtt_sim, 1, score_fn, na.rm=TRUE)
        score_sim[is.infinite(score_sim)] <- NA
        Fn <- stats::ecdf(c(score_sim, Inf, -Inf))
        if(d=="up"){
          pval <- 1 - Fn(score)
        }else if(d=="down"){
          pval <- Fn(score)
        }
      }else{
        stop("Method must be 'exact', 'asymptotic', or 'monte_carlo'")
      }
      fdr <- stats::p.adjust(pval, method=adjust.method)
      res_tmp <- data.frame(score=score, p=pval, FDR=fdr)
      colnames(res_tmp) <- paste(groups[i], d, colnames(res_tmp), sep=".")
      res[[i]] <- res_tmp
    }
    res <- Reduce(cbind, res)
    res <- res[rownames(dds), ]
    res_all[[d]] <- res
  }

  res_all <- Reduce(cbind, res_all)
  if(add.means){
    mat <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
    mat_avg <- sapply(groups, function(g) rowMeans(mat[, grp==g, drop=FALSE]))
    colnames(mat_avg) <- paste0(groups, ".avg")
    res_all <- cbind(mat_avg[rownames(res_all), ], res_all)
  }
  return(res_all)
}
