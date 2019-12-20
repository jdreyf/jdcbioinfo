#' Parse kallisto log files
#'
#' Parse kallisto log files.
#'
#' @param name  the file name of result.
#' @param path  the full path name for the kallisto log files.
#' @return Data frame with summary of kallisto alignment.
#' @export

parse_kallisto_logs <- function(name="kallisto_summary", path="."){

  out_files <- dir(path=path, pattern="\\.[0-9]+\\.out$")

  tab <- list()
  for (f in out_files) {
    samp <- gsub("\\.[0-9]+\\.out$", "", f)
    lines <- readLines(f)
    res <- grep("\\[quant\\] processed", lines, value=TRUE)
    res <- unlist(strsplit(res, split=" "))[c(3, 5)]
    res <- as.numeric(gsub(",", "", res))
    res[3] <- round(res[2]/res[1] * 100, 2)

    len <- grep("estimated average fragment length", lines, value=TRUE)
    if(length(len)==0) {
      len <- NA
    } else {
      len <- as.numeric(trimws(unlist(strsplit(len, split=":"))[2]))
    }
    res[4] <- len

    names(res) <- c("Total reads", "Pseudoaligned reads", "Percent of aligned", "Estimated average fragment length")
    tab[[samp]] <- res
  }

  tab <- sapply(tab, FUN=function(v) v)
  if (!is.na(name)) utils::write.csv(tab, paste0(name, ".csv"))

  return(invisible(tab))
}
