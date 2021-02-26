#' Impute missing values by random numbers drawn from a normal distribution
#'
#' Impute missing values by random numbers drawn from a normal distribution that has a down-shifted mean
#' and shrunken standard deviation from the sample distribution. This is meant to be similar to imputation
#' in the Perseus software.
#'
#' @param object Data frame or matrix containing filtered and log-transformed data.
#' @param width Scale factor for the standard deviation of imputed distribution relative to the sample standard deviation.
#' @param downshift Down-shifted the mean of imputed distribution from the sample mean, in units of sample standard deviation.
#' @param seed Random seed
#' @return A matrix with imputed values
#' @export

impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {

  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")

  set.seed(seed)
  object <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
    })
  return(object)
}
