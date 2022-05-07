##' Enrichment analyzer using \code{\link[limma]{roast}} with function \code{fry}
##'
##'
##' @title multi_enrichFLY
##' @param contrast.v Named vector of contrasts, passed to \code{\link[limma]{makeContrasts}}.
##' @inheritParams enrichFLY
##' @return  A list of 3 \code{compareClusterResult} instance for up-regulated, down-regulated, and mixed-regulated enriched gene sets.
##' @importClassesFrom DOSE compareClusterResult
##' @export

multi_enrichFLY <- function(object, G, annot, sep.str = " /// ", symbol.col = "symbol", grp = NULL, contrast.v = ncol(design), design = NULL,
                      weights = NULL, trend = FALSE, block = NULL, correlation = NULL, adjust.method = c("BH", "none"),
                      min.nfeats = 3, max.nfeats = 1000, pvalueCutoff = 0.25, qvalueCutoff = 1) {
  enrichResultList <- lapply(contrast.v, function(contrast) {
    enrichFLY(object, G = G, annot = annot, sep.str = sep.str, symbol.col = symbol.col, grp = grp, contrast = contrast, design = design,
              weights = weights, trend = trend, block = block, correlation = correlation, adjust.method = adjust.method,
              min.nfeats = min.nfeats, max.nfeats = max.nfeats, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)})

  upEnrichResultList <- lapply(enrichResultList, "[[", "Up")
  downEnrichResultList <- lapply(enrichResultList, "[[", "Down")
  mixEnrichResultList <- lapply(enrichResultList, "[[", "Mixed")

  up <- clusterProfiler::merge_result(upEnrichResultList)
  down <- clusterProfiler::merge_result(downEnrichResultList)
  mix <- clusterProfiler::merge_result(mixEnrichResultList)

  return(list(Up = up, Down = down, Mixed = mix) )
}

