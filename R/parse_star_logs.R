#' Parse star log files
#'
#' Parse star log files.
#'
#' @inheritParams parse_kallisto_logs
#' @return Data frame with summary of star alignment.
#' @export
#'
parse_star_logs <- function(name="star_summary", path=".") {

  out_files <- dir(path=path, pattern="Log\\.final\\.out$", recursive=TRUE, full.names=TRUE)
  stopifnot(length(out_files) > 1)
  tab <- list()
  for (f in out_files) {
    samp <- gsub(".*/", "", dirname(f))
    res <- utils::read.delim(f, header=FALSE)
    res[,1] <- trimws(gsub("\\|", "", res[,1]))
    res <- stats::setNames(trimws(res[,2]), nm=res[,1])
    tab[[samp]] <- res
  }
  tab <- sapply(tab, FUN=function(v) v)
  if (!is.na(name)) utils::write.csv(tab, paste0(name, ".csv"))

  return(invisible(tab))
}
