#' Make file links forR Markdown
#'
#'Make file links forR Markdown
#'
#' @inheritParams base::dir
#' @return Vector of file links.
#' @export

make_file_links <- function(path=".", pattern=NULL, recursive=TRUE, full.names=TRUE) {
  full.path <- dir(path=path, pattern=pattern, recursive=recursive, full.names=full.names)
  file.nm <- basename(full.path)
  links <- paste0("[", file.nm, "](", full.path, ")")
  return(links)
}
