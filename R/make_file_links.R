#' Make file links forR Markdown
#'
#'Make file links forR Markdown
#'
#' @inheritParams base::dir
#' @return Vector of file links.
#' @export

make_file_links <- function(path=".", pattern=NULL, recursive=TRUE, full.names=TRUE) {
  full.path <- dir(path=path, pattern=pattern, recursive=recursive, full.names=full.names)

  file.nm <- full.path
  idxs <- gregexpr("/", file.nm)
  mx <- max(sapply(idxs, length))

  for(i in seq_len(mx)){
    file.nm.tmp <- sub("[^/]*/", "", file.nm)
    if (sum(duplicated(file.nm.tmp)) > 0) break
    file.nm <- file.nm.tmp
  }

  links <- paste0("[", file.nm, "](", full.path, ")")
  return(links)
}
