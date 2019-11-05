#' Filter features
#'
#' Filter features that have more than a proportion missing values in one or all groups.
#'
#' @param object A dat.frame or matrix containing data for filtering.
#' @param grp A character vector dictating the grouping.
#' @param min.prop  The propotion of valid values for each group for retention.
#' @param at.least.one TRUE means to keep the row if min.prop is met for at least one group,
#' FALSE means min.prop must be met across all groups for retention.
#' @param only.keep.vector TRUE for returning only the logical vector indicates the rows to keep,
#' FALSE for returning the matix after filtering.
#' @return Logical vector or data matirx.
#' @export

filter_valids <- function(object, grp, min.prop=2/3, at.least.one=FALSE, only.keep.vector=FALSE){

  stopifnot(ncol(object)==length(grp), colnames(object)==names(grp))

  if (!is.matrix(object)) object <- as.matrix(object)

  grp_nms <- unique(grp)
  grp_filter <- vapply(grp_nms, FUN=function(nm) {
    object_tmp <- object[, grp==nm]  # Extract columns of interest
    sums <- rowSums(is.finite(object_tmp)) # Count the number of valid values for each group
    sums >=round(min.prop * ncol(object_tmp))   # Calculates whether min.prop requirement is met
    }, FUN.VALUE=logical(nrow(object)))

  if (at.least.one) {
    filter <- apply(grp_filter, 1, any)
  } else {
    filter <- apply(grp_filter, 1, all)
  }
  if (only.keep.vector) return(filter)

  object <- object[filter, ]
  return(object)
}
