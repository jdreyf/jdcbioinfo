#' Add beta-values to table with means of M-values
#'
#' Add beta-values to table with means of M-values, where beta-values are proportion methylated and M-values are their logit.
#'
#' @param tt Data frame with columns that have means of M-values and end with `.avg`.
#' @return Data frame with additional columns of mean beta-values.
#' @export

add_avg_beta <- function(tt){
  stopifnot(any(grepl(pattern = ".avg", x=colnames(tt), fixed = TRUE)), !grepl(pattern = "_beta.avg", x=colnames(tt), fixed = TRUE))
  tt <- tt |> dplyr::select(tidyselect::ends_with("avg")) |>
    dplyr::mutate(dplyr::across(.cols = tidyselect::everything(), .fn = boot::inv.logit)) |>
    dplyr::rename_with(.fn = \(x) sub(".avg", "_beta.avg", x, fixed = TRUE)) |>
    dplyr::bind_cols(tt)
  tt
}
