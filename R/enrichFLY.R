##' Enrichment analyzer using \code{\link[limma]{roast}} with function \code{fry}
##'
##'
##' @title enrichFLY
##' @param object Matrix-like data object containing log-ratios or log-expression values, with rows corresponding to
##' features (e.g. genes) and columns to samples. Must have row names that are non-duplicated and non-empty.
##' @param G Gene set list as returned from \code{\link{read_gmt}}.
##' @param annot annotation for the features that has a column of the same type as in gene set list `G`.
##' @param sep.str strings that separates symbols if there are multiple symbols for a feature.
##' @param symbol.col column name or index for the symbol column in `annot`.
##' @param grp Vector of sample groups. These must be valid variable names in R and the same length as
##' \code{ncol(object)}.
##' @param contrast contrast for which the test is required. Can be an integer specifying a column of `design`,
##' or the name of a column of `design`, or a numeric contrast vector of length equal to the number of columns of `design`.
##' @param design Design matrix of the experiment, with rows corresponding to samples and columns to coefficients to be
##' estimated.
##' @param weights Non-negative observation weights. Can be a numeric matrix of individual weights of same size as the
##' \code{object}, or a numeric vector of sample weights with length \code{ncol(object)}, or a numeric vector of gene
##' weights with length equal to \code{nrow(object)}.
##' @param trend Logical; should an intensity-trend be allowed for the prior variance? Default is that the prior variance
##' is constant.
##' @param block Vector specifying a blocking variable on the samples. Has length = \code{ncol(object)}.
##' Must be \code{NULL} if \code{ndups > 1}.
##' @param correlation Inter-duplicate or inter-technical replicate correlation. Must be given if
##' \code{ndups>1} or \code{!is.null(block)}.
##' @param adjust.method Method used to adjust the p-values for multiple testing. Options,
##' @param min.nfeats Minimum number of features (e.g. genes) needed in a gene set for testing.
##' @param max.nfeats Maximum number of features (e.g. genes) needed in a gene set for testing.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param qvalueCutoff Cutoff value of qvalue.
##' @return  A list of 3 \code{enrichResult} instance for up-regulated, down-regulated, and mixed-regulated enriched gene sets.
##' @importClassesFrom DOSE enrichResult
##' @export

enrichFLY <- function(object, G, annot, sep.str = " /// ", symbol.col = "symbol", grp = NULL, contrast = ncol(design), design = NULL,
                      weights = NULL, trend = FALSE, block = NULL, correlation = NULL, adjust.method = c("BH", "none"),
                      min.nfeats = 3, max.nfeats = 1000, pvalueCutoff = 0.25, qvalueCutoff = 1) {
  stopifnot(!is.null(dim(object)), !is.null(rownames(object)),
            !is.null(colnames(object)),
            ncol(object) > 1,
            !is.null(design) || !is.null(grp),
            length(weights) != 1 || is.null(weights),
            length(weights) <=  1 || (is.numeric(weights) && all(weights >= 0) && !all(is.null(weights))),
            length(weights) <= 1 || all(dim(weights) ==  dim(object)) || length(weights) == nrow(object) || length(weights) == ncol(object),
            rownames(object) %in% rownames(annot))

  if (!is.null(block) && is.null(correlation)) stop("!is.null(block), so correlation must not be NULL.")

  adjust.method <- match.arg(adjust.method, choices = c("BH", "none"))
  index <- ezlimma:::g_index(G = G, object = object, min.nfeats = min.nfeats,  max.nfeats = max.nfeats)

  if (is.null(design)) {
    stopifnot(ncol(object) == length(grp))
    design <- stats::model.matrix(~0 + grp)
    colnames(design) <- sub("grp", "", colnames(design), fixed = TRUE)
  }
  contr.mat <- limma::makeContrasts(contrasts = contrast, levels = design)
  colnames(contr.mat) <- names(contrast)

  if (!is.matrix(object)) {
    if (!is.null(object$weights)) {
      if (!is.null(weights)) {
        warning("object$weights are being ignored")
      } else {
        if (is.null(dimnames(object$weights))) dimnames(object$weights) <- dimnames(object)
        weights <- object$weights
      }
    }
  }

  #	Covariate for trended eBayes
  covariate <- NULL
  if(trend) covariate <- rowMeans(as.matrix(object))

  #	Compute effects matrix (with df.residual+1 columns)
  Effects <- limma:::.lmEffects(object, design = design, contrast = contr.mat, weights = weights, block = block, correlation = correlation)

  # Divide out genewise standard deviations
  # Estimate genewise sds robustly
  OK <- requireNamespace("statmod", quietly = TRUE)
  if(!OK) stop("statmod package required but isn't installed (or can't be loaded)")
  gq <- statmod::gauss.quad.prob(128,"uniform")
  df.residual <- ncol(Effects) - 1
  Eu2max <- sum((df.residual + 1) * gq$nodes^df.residual * stats::qchisq(gq$nodes, df = 1) * gq$weights)
  u2max <- apply(Effects^2, 1, max)
  s2.robust <- (rowSums(Effects^2) - u2max) / (df.residual + 1 - Eu2max)

  # Estimate hyperparameters from residual variances but apply squeezing to robust variances
  s2 <- rowMeans(Effects[, -1, drop = FALSE]^2)
  fit <- limma::fitFDist(s2, df1 = df.residual, covariate = covariate)
  df.prior <- fit$df2

  s2.robust <- limma:::.squeezeVar(s2.robust, df = 0.92 * df.residual, var.prior = fit$scale, df.prior = df.prior)
  Effects <- Effects/sqrt(s2.robust)

  # Gene ID
  geneid <- rownames(Effects)

  # Fry effect
  nfeatrues <- nrow(Effects)
  neffects <- ncol(Effects)
  df.residual <- neffects - 1L

  #	Check index
  if(is.null(index)) index <- list(set1 = 1L:nfeatrues)
  if(!is.list(index)) index <- list(set1=index)
  nsets <- length(index)
  if(nsets == 0L) stop("index is empty")
  if(is.null(names(index))) {
    names(index) <- paste0("set",formatC(1L:nsets, width = 1L + floor(log10(nsets)), flag = "0"))
  } else {
     if(anyDuplicated(names(index))) stop("Gene sets don't have unique names", call. = FALSE)
  }

  # Overall statistics
  MeanEffects <- colMeans(Effects)

  t.cutoff <- stats::qt(0.025, df = df.residual, lower.tail = FALSE)
  t.stat.all <- sapply(Effects[,1], function(x) x / sqrt(mean(MeanEffects[-1]^2)))
  sig.gene.all <- rep_len(0L, length.out = length(t.stat.all))
  names(sig.gene.all) <- names(t.stat.all)
  sig.gene.all[t.stat.all > t.cutoff] <- 1L
  sig.gene.all[t.stat.all < -t.cutoff] <- -1L

  # Set statistics
  NGenes <- rep_len(0L, length.out = nsets)
  PValue.Mixed <- t.stat <- rep_len(0, length.out =  nsets)
  sig.gene.set <- list()
  for (i in 1:nsets) {
    iset <- index[[i]]
    if(is.factor(iset)) iset <- as.character(iset)
    if(is.character(iset)) iset <- which(geneid %in% iset)
    EffectsSet <- Effects[iset,,drop=FALSE]

    MeanEffectsSet <- colMeans(EffectsSet)

    t.stat.set <- sapply(EffectsSet[,1], function(x) x / sqrt(mean(MeanEffectsSet[-1]^2)))
    sig.gene.set[[i]] <- rep_len(0L, length.out = length(t.stat.set))
    names(sig.gene.set[[i]]) <- names(t.stat.set)
    sig.gene.set[[i]][t.stat.set > t.cutoff] <- 1L
    sig.gene.set[[i]][t.stat.set < -t.cutoff] <- -1L

    t.stat[i] <- MeanEffectsSet[1] / sqrt(mean(MeanEffectsSet[-1]^2))
    NGenes[i] <- nrow(EffectsSet)

    if(NGenes[i] > 1) {
      SVD <- svd(EffectsSet, nu = 0)
      A <- SVD$d^2
      d1 <- length(A)
      d <- d1 - 1L
      beta.mean <- 1/d1
      beta.var <- d/d1/d1/(d1/2 + 1)
      Fobs <- (sum(EffectsSet[,1]^2) - A[d1]) / (A[1] - A[d1])
      Frb.mean <- (sum(A) * beta.mean - A[d1]) / (A[1] - A[d1])
      COV <- matrix(-beta.var/d, d1, d1)
      diag(COV) <- beta.var
      Frb.var <- (A %*% COV %*% A ) / (A[1] - A[d1])^2
      alphaplusbeta <- Frb.mean*(1 - Frb.mean)/Frb.var - 1
      alpha <- alphaplusbeta*Frb.mean
      beta <- alphaplusbeta - alpha
      PValue.Mixed[i] <- stats::pbeta(Fobs, shape1 = alpha, shape2 = beta,lower.tail = FALSE)
    }
  }

  Direction <- rep_len("Up", length.out = nsets)
  Direction[t.stat < 0] <- "Down"
  PValue <- 2*stats::pt(-abs(t.stat),df = df.residual)
  PValue.Mixed[NGenes == 1] <- PValue[NGenes == 1]

  #	FDR
  if(nsets > 1 && adjust.method == "BH") {
    FDR <- stats::p.adjust(PValue, method = "BH")
    FDR.Mixed <- stats::p.adjust(PValue.Mixed, method = "BH")
  } else {
    FDR <- PValue
    FDR.Mixed <- PValue.Mixed
  }

  # qval
  qobj <- tryCatch(qvalue::qvalue(p = PValue, lambda = 0.05, pi0.method = "bootstrap"), error=function(e) NULL)
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }

  qobj.mix <- tryCatch(qvalue::qvalue(p = PValue.Mixed, lambda = 0.05, pi0.method = "bootstrap"), error=function(e) NULL)
  if (class(qobj.mix) == "qvalue") {
    qvalues.mix <- qobj.mix$qvalues
  } else {
    qvalues.mix <- NA
  }

  # significant genes up/down
  sig.gene.reg <- lapply(seq_along(sig.gene.set), function(i) {
    if (Direction[i] == "Up") {
      names(sig.gene.set[[i]])[sig.gene.set[[i]] == 1]
    } else {
      names(sig.gene.set[[i]])[sig.gene.set[[i]] == -1]
    }
  })
  count.reg  <- sapply(sig.gene.reg, length)
  sig.gene.reg <- lapply(sig.gene.reg, function(ids) annot[ids, symbol.col])
  sig.gene.reg <- lapply(sig.gene.reg, function(genes) {
   unique(unlist(strsplit(genes, split = sep.str)))
  })
  sig.gene.reg <- sapply(sig.gene.reg, paste, collapse = "/")

  # significant genes mixed
  sig.gene.mix <- lapply(seq_along(sig.gene.set), function(i) {
      names(sig.gene.set[[i]])[sig.gene.set[[i]] != 0]
  })
  count.mix  <- sapply(sig.gene.mix, length)
  sig.gene.mix <- lapply(sig.gene.mix, function(ids) annot[ids, symbol.col])
  sig.gene.mix <- lapply(sig.gene.mix, function(genes) {
    unique(unlist(strsplit(genes, split = sep.str)))
  })
  sig.gene.mix <- sapply(sig.gene.mix, paste, collapse = "/")

  tab.reg <- data.frame(ID          = names(index),
                        Description = names(index),
                        GeneRatio   = paste0(count.reg, "/", NGenes),
                        BgRatio     = ifelse(Direction == "Up", paste0(sum(sig.gene.all == 1), "/", length(sig.gene.all)),
                                             paste0(sum(sig.gene.all == -1), "/", length(sig.gene.all))),
                        pvalue      = PValue,
                        p.adjust    = FDR,
                        qvalue      = qvalues,
                        geneID      = sig.gene.reg,
                        Count       = count.reg)

  tab.up <- tab.reg[Direction == "Up", ]
  tab.down <- tab.reg[Direction == "Down", ]
  tab.up <- tab.up[order(tab.up$pvalue), ]
  tab.down <- tab.down[order(tab.down$pvalue), ]

  tab.mix <- data.frame(ID          = names(index),
                        Description = names(index),
                        GeneRatio   = paste0(count.mix, "/", NGenes),
                        BgRatio     = paste0(sum(sig.gene.all != 0), "/", length(sig.gene.all)),
                        pvalue      = PValue.Mixed,
                        p.adjust    = FDR.Mixed,
                        qvalue      = qvalues.mix,
                        geneID      = sig.gene.mix,
                        Count       = count.mix)
  tab.mix <- tab.mix[order(PValue.Mixed), ]

  up <- methods::new("enrichResult",
                      result         = tab.up,
                      pvalueCutoff   = pvalueCutoff,
                      pAdjustMethod  = adjust.method,
                      qvalueCutoff   = qvalueCutoff,
                      gene           = names(sig.gene.all)[sig.gene.all == 1],
                      universe       = names(sig.gene.all),
                      geneSets       = index,
                      organism       = "UNKNOWN",
                      keytype        = "UNKNOWN",
                      ontology       = "UNKNOWN",
                      readable       = FALSE)

  down <- methods::new("enrichResult",
                        result         = tab.down,
                        pvalueCutoff   = pvalueCutoff,
                        pAdjustMethod  = adjust.method,
                        qvalueCutoff   = qvalueCutoff,
                        gene           = names(sig.gene.all)[sig.gene.all == -1],
                        universe       = names(sig.gene.all),
                        geneSets       = index,
                        organism       = "UNKNOWN",
                        keytype        = "UNKNOWN",
                        ontology       = "UNKNOWN",
                        readable       = FALSE)

  mix <- methods::new("enrichResult",
                       result         = tab.mix,
                       pvalueCutoff   = pvalueCutoff,
                       pAdjustMethod  = adjust.method,
                       qvalueCutoff   = qvalueCutoff,
                       gene           = names(sig.gene.all)[sig.gene.all != 0],
                       universe       = names(sig.gene.all),
                       geneSets       = index,
                       organism       = "UNKNOWN",
                       keytype        = "UNKNOWN",
                       ontology       = "UNKNOWN",
                       readable       = FALSE)

  # subset
  up <- DOSE:::get_enriched(up)
  if (nrow(up@result) == 0) warning("No up-regualted gene sets are enriched")
  down <- DOSE:::get_enriched(down)
  if (nrow(down@result) == 0) warning("No down-regualted gene sets are enriched")
  mix <- DOSE:::get_enriched(mix)
  if (nrow(mix@result) == 0) warning("No mixed-regualted gene sets are enriched")

  return(list(Up = up, Down = down, Mixed = mix))
}

