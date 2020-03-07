#' Rand seq
#'
#' Rand seq.
#'
#' @param BSgenome BSgenome object.
#' @param nsamp Numeric number of samples.
#' @param flank Numeric length of flank on each side of region.
#' @inheritParams ezlimma::roast_contrasts

sample_bsgenome_regions <- function(BSgenome, nsamp, flank, seed = 123){
  bed.samp <- vector("list", nsamp)
  set.seed(seed)
  for(samp.ind in seq_along(bed.samp)){
    prob <- BSgenome@seqinfo@seqlengths / sum(as.numeric(BSgenome@seqinfo@seqlengths))
    idx <- sample(length(BSgenome@seqinfo), 1, prob = prob)
    chr <- BSgenome@seqinfo@seqnames[idx]
    seqlength <- BSgenome@seqinfo@seqlengths[idx]
    start <- sample(seqlength - 2*flank, 1)
    end <- start + 2*flank
    name <- paste0("random_seq_", samp.ind)
    bed.samp[[samp.ind]] <- data.frame(chr=chr, start = start, end = end, name = name)
  }
  bed.samp <- Reduce(rbind, bed.samp)
  rownames(bed.samp) <- bed.samp$name
  gr.samp <- methods::as(bed.samp, "GRanges")

  # test
  seqs.samp <- BSgenome::getSeq(BSgenome, gr.samp)
  stopifnot(names(seqs.samp) == names(gr.samp))

  return(seqs.samp)
}
