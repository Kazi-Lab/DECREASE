#' Create Bliss independence matrix from dose-response matrix
#'
#' @details Values are limited between 0 and 100 before calculations are made
#'
#' @param mat A dose response matrix, where values are % inhibition
#' @returns A Bliss independence matrix, where values are % inhibition
#' @examples rnorm(25, mean = 50, sd = 50, 5, 5) |> calculate_bliss()
#' @export
calculate_bliss <- function(mat) {
  # Limit between 0 and 100
  mat <- matrix(pmax(0, pmin(mat, 100)), nrow(mat), ncol(mat))

  # Convert from inhibition to survival probability
  surv <- (1 - mat / 100)

  # Take outer product of first row and column
  bliss <- surv[2:nrow(surv), 1] %o% surv[1, 2:ncol(surv)]

  # Replace combination responses with Bliss independece reference values
  surv[2:nrow(surv), 2:ncol(surv)] <- bliss

  # Convert back to inhibition
  (1 - surv) * 100
}
