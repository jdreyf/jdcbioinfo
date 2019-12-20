#' Parse cutadapt log files
#'
#' Parse cutadapt log files.
#'
#' @inheritParams parse_kallisto_logs
#' @return Data frame with summary of cutadapt.
#' @export
#'
parse_cutadapt_logs <- function(name="cutadapt_summary", path=".") {

  out_files <- dir(path=path, pattern="\\.[0-9]+\\.out$")
  stopifnot(length(out_files) > 1)
  res <- list()
  for (f in out_files) {
    nm <- gsub("\\.[0-9]+\\.out", "", basename(f))
    lines <- readLines(f)
    idx1 <- grep("Total read.*processed", lines)
    idx2 <- grep("written \\(passing filters\\)", lines)
    lines <- lines[idx1:idx2]
    lines <- strsplit(lines, split = ":")
    lines <- lapply(lines, FUN=function(line) gsub(",|\\(|\\)|%", "", trimws(line)))
    lines <- lapply(lines, function(line) c(line[1], unlist(strsplit(line[2], split = " +"))))
    dfs <- lapply(lines, FUN=function(line) {
      if (length(line) == 2) {
        df <- data.frame(rn = line[1], stat = as.numeric(line[2]), stringsAsFactors = FALSE)
        } else {
          stopifnot(length(line) == 3)
          df <- data.frame(rn = line[1], stat = as.numeric(line[2:3]), stringsAsFactors = FALSE)
          df$rn[2] <- paste(df$rn[2], "(%)")
        }
        colnames(df)[2] <- nm
        df })
        res[[nm]] <- Reduce(rbind, dfs)
    }
    res <- Reduce(function(x,y){merge(x , y, all = TRUE, sort = FALSE)}, res)
    res <- data.frame(res, row.names = "rn")

    if (!is.na(name)) utils::write.csv(res, paste0(name, ".csv"))

    return(invisible(res))
}
