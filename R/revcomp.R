#' converts a DNA sequence into its reverse-complement counterpart
#'
#' converts a DNA sequence into its reverse-complement counterpart.
#'
#' @param x a character vector of strings of DNA sequences
#' @return a character vector of strings of the reverse-complement DNA sequences
#' @export

revcomp <- function(x){
  sapply(lapply(strsplit(chartr("ACTGactg", "TGACtgac", x), NULL), rev), paste, collapse="")
}
