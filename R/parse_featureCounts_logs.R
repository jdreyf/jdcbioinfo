#' Parse featureCounts log file
#'
#' Parse featureCounts log file.
#'
#' @param file  the log file name including the full path.
#' @param name  the file name of result.
#' @return Data frame with summary of featureCounts.
#' @export

parse_featureCounts_logs <- function(file, name="featureCounts_summary"){

  lines <- readLines(file)

  samps <- grep("Process BAM file", lines, value = TRUE)
  samps <- gsub(".*( |/)", "", gsub("/Aligned\\..*out\\.bam.*", "", samps))

  tot.num <- grep("Total fragments", lines, value = TRUE)
  tot.num <- sapply(strsplit(trimws(gsub("\\|\\|", "", tot.num)), split = " : "), FUN="[", 2)

  mapped <- grep("Successfully assigned fragments", lines, value = TRUE)
  mapped <- strsplit(trimws(gsub("\\|\\|", "", mapped)), split = " ")
  mapped.num <- sapply(mapped, FUN="[", 5)
  mapped.percent <- gsub("\\(|\\)", "", sapply(mapped, FUN="[", 6))

  tab <- t(cbind(tot.num, mapped.num, mapped.percent))
  dimnames(tab) <- list(c("Total fragments", "Successfully assigned fragments", "Percentage of assigned fragments"),
                         samps)

  if (!is.na(name)) utils::write.csv(tab, paste0(name, ".csv"))

  return(invisible(tab))
}
