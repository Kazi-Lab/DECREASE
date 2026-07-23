#' Mark TRUE if values are outside of range
#'
#' An outlier here is defined as any value that is outside of the (25%, 75%)
#' quartile range, plus (or minus) iqrs*IQR. If iqrs = 0, values outside of
#' (25%, 75%) will be marked as outliers. If iqrs = 1, values outside of ((iqrs
#' * IQR) - 25%), (iqrs * IQR) + 75%) will be marked as outliers
#'
#' @param values Values to be marked as outliers
#' @param iqrs The number of additional IQRS to extend the range by. Can be
#'   fractional.
#' @returns A vector of booleans
#' @examples
#' set.seed(1)
#' values <- rnorm(100, sd = 10)
#'
#' mark_outliers(values, iqrs = 0)
#' mark_outliers(values, iqrs = 1.5)
#' @export
mark_outliers <- function(values, iqrs) {
  quantiles <- stats::quantile(values, probs = c(0.25, 0.75), na.rm = TRUE)
  additional_range <- iqrs * stats::IQR(values, na.rm = TRUE)
  ll <- quantiles[1] - additional_range
  ul <- quantiles[2] + additional_range
  values < ll | values > ul
}

detect_outliers <- function(mat, fcurve, fit_1, fit_2) {
  d1 <- fit_1$response
  d2 <- fit_2$response

  bliss_mat <- calculate_bliss(mat)

  if (fcurve) {
    mat[, 1] <- d1
    mat[1, ] <- d2
  } else {
    mat[, 1] <- (mat[, 1] + d1) / 2
    mat[1, ] <- (mat[1, ] + d2) / 2
  }

  dev_d1 <- abs(d1 - (mat[, 1]))
  dev_d2 <- abs(d2 - (mat[1, ]))
  dev <- abs(bliss_mat - mat)
  outliers <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

  # check deviations with Bliss
  outliers[-1, -1] <- mark_outliers(
    abs(dev[-1, -1] - stats::median(dev[-1, -1], na.rm = TRUE)),
    5
  ) &
    (dev[-1, -1] > 25)

  outliers[, 1] <- (apply(outliers, 2, any, na.rm = TRUE) & dev_d1 > 10) |
    dev_d1 > 15
  outliers[1, ] <- (apply(outliers, 1, any, na.rm = TRUE) & dev_d2 > 10) |
    dev_d1 > 15
  # BUG: dev_d1 > 15 is used for both statements. This is likely a mistake.

  outliers
}
