#' Fit dose response data
#'
#' @details Broadly, this function:
#' 1. Converts dose into log10 dose
#' 2. Attempts a few preliminary fits to get good starting parameters
#' 3. Checks if parameters seem unrealistic and nudges them accordingly
#' 4. Does a final fit with these parameters
#'
#' @param data data.frame. Has two columns:
#' * dose: numeric, contains untransformed doses
#' * response: numeric, contains responses as inhibition
#' @export
fit_dose_response <- function(data) {
  data$log_dose <- log10(data$dose)
  data <- data[order(data$dose), ]

  data <- break_response_ties(data)

  estimated_parameters <- estimate_parameters(data)

  coef_estim <- stats::coef(estimated_parameters)
  names(coef_estim) <- c("SLOPE", "MIN", "MAX", "IC50")

  cleaned_parameters <- clean_parameters(coef_estim, data)

  nls_fit(data, cleaned_parameters)
}

#' Break response ties
#'
#' If there are any ties, breaks them by adding an escalating 0.01 to EVERY
#' response
#'
#' @inheritParams fit_dose_response
#'
#' @returns input with responses broken, if any
#'
#' @examples
#' no_ties <- data.frame(dose = 1:3, response = 1:3)
#' break_response_ties(no_ties)
#'
#' ties <- data.frame(dose = 1:3, response = c(1, 1, 3))
#' break_response_ties(ties)
break_response_ties <- function(data) {
  if (any(duplicated(data$response))) {
    data$response <- seq(
      from = 0,
      length.out = length(data$response),
      by = 0.01
    ) +
      data$response
  }
  data
}

drc_fit_with_fn <- function(data, fn, names, errorm = TRUE) {
  drc::drm(
    response ~ log_dose,
    data = data,
    fct = fn(names = names),
    logDose = 10,
    control = drc::drmc(errorm = errorm)
  )
}

#' Try to fit data to get initial parameters
#'
#' Will attempt with log-logistic drc::LL.4. If fail or warn, falls back to
#' drc::L.4.
#'
#' @inheritParams fit_dose_response
#'
#' @returns a object of class `drc`
estimate_parameters <- function(x) {
  names <- c("SLOPE", "MIN", "MAX", "IC50")
  tryCatch(
    drc_fit_with_fn(x, drc::LL.4, names, errorm = FALSE),
    warning = \(w) drc_fit_with_fn(x, drc::L.4, names),
    error = \(e) drc_fit_with_fn(x, drc::L.4, names)
  )
}

#' Detect and nudge suspicious looking parameters
#'
#' @inheritParams fit_dose_response
#'
#' @param p A named vector of coefficients of a fit, with names "SLOPE", "MIN",
#'   "MAX", and "IC50"
clean_parameters <- function(p, data) {
  # p = parameters
  #
  # I know it's frowned upon to abbreviate variables, but less text makes it less
  # scary to look at
  max_dose <- max(data$dose, na.rm = TRUE)
  min_dose <- min(data$dose, na.rm = TRUE)
  max_response <- max(data$response, na.rm = TRUE)
  min_response <- min(data$response, na.rm = TRUE)
  max_log_dose <- max(data$log_dose, na.rm = TRUE)
  p <- as.list(p)

  # drc's LL.4/L.4 curves are parameterized atypically where an increasing
  # dose-response curve yields a negative slope.
  p$SLOPE <- -p$SLOPE

  if (p$MAX <= p$MIN || p$IC50 > max_dose) {
    p$IC50 <- max_dose
  }

  # Previously this logic would check if it was below 0 and set it to the min
  # dose, then check if it was below 0 and set it to the mean dose. This would
  # only fire if the min dose was negative, which should never happen
  if (p$IC50 < 0) {
    p$IC50 <- min_dose
  }

  # This will likely cause problems if the previous statement fired, because the
  # min_dose is likely to be 0. However, this is what was previously written.
  p$IC50 <- log10(p$IC50)

  # NOTE: OG version didn't use na.rm
  if (p$IC50 < min(data$log_dose)) {
    p$IC50 <- max(data$log_dose)
  }

  if (all(data$response < 0)) {
    p$IC50 <- max_log_dose
  }

  p$MIN <- 0
  p$MAX <- min(max_response, 100)

  min_lower <- min(99, max(min_response, 0))
  if (p$MAX < 0) {
    p$MAX <- 100
  }
  max_lower <- ifelse(
    max_response <= 100 && max_response >= 0,
    max_response,
    100
  )

  run_avg <- caTools::runmean(data$response, 10)
  max_run_avg_index <- which.max(run_avg)
  max_upper <- ifelse(
    max_run_avg_index != length(run_avg),
    run_avg[max_run_avg_index],
    p$MAX
  )

  if (any(data$response > max_upper)) {
    max_upper <- mean(data$response[data$response > max_upper]) + 5
  }
  if (max_upper < 0) {
    max_upper <- p$MAX
  }
  max_upper <- min(max_upper, 100)
  if (max_lower > max_upper) {
    max_upper <- p$MAX
  }

  mean_resp_last <- mean(utils::tail(data$response, 2), na.rm = TRUE)
  if (mean_resp_last < 25) {
    p$IC50 <- max_log_dose
  } else if (mean_resp_last < 60) {
    p$IC50 <- mean(data$log_dose, na.rm = TRUE)
  }
  if (mean(data$response[1:3], na.rm = TRUE) < 5) {
    p$IC50 <- max_log_dose
  }
  if (p$MIN == p$MAX) {
    p$MAX <- p$MAX + 0.001
  }

  list(
    params = unlist(p),
    max_upper = max_upper,
    max_lower = max_lower,
    min_lower = min_lower,
    mean_resp_last = mean_resp_last
  )
}

#' Fit dose response data with best-guess starting parameters
#'
#' @inheritParams fit_dose_response
#'
#' @param cleaned_parameters list containing:
#' * params: a named vector of fit coefficients
#' * max_lower, min lower, and max_upper: additional constrains to provide nls
#' if slope is shallow
#' @returns nls fit
nls_fit <- function(data, cleaned_parameters) {
  coef_estim <- cleaned_parameters$params
  max_lower <- cleaned_parameters$max_lower
  min_lower <- cleaned_parameters$min_lower
  max_upper <- cleaned_parameters$max_upper

  form <- response ~ MIN + (MAX - MIN) / (1 + (10^(SLOPE * (IC50 - log_dose))))

  start <- c(
    SLOPE = 1,
    MIN = coef_estim["MIN"][[1]],
    MAX = coef_estim["MAX"][[1]],
    IC50 = coef_estim["IC50"][[1]]
  )

  ic50_median_dose_start <- start
  ic50_median_dose_start[["IC50"]] <- stats::median(data$log_dose)

  lower <- c(SLOPE = 0, MIN = 0, MAX = max_lower, IC50 = min(data$log_dose))
  upper <- c(SLOPE = 4, MIN = 0, MAX = 100, IC50 = max(data$log_dose))

  # Settings for nls
  control <- list(warnOnly = TRUE, minFactor = 1 / 2048)
  fit_1 <- tryCatch(
    {
      stats::nls(
        form,
        data,
        start,
        control,
        "port",
        lower = lower,
        upper = upper
      )
    },
    error = \(e) {
      minpack.lm::nlsLM(form, data, start, lower = lower, upper = upper)
    }
  )

  fit_2 <- tryCatch({
    stats::nls(
      form,
      data,
      ic50_median_dose_start,
      control,
      "port",
      lower = lower,
      upper = upper
    )
  })

  best_fit <- choose_best_fit(fit_1, fit_2)

  # In the case of a shallow slope, fit again with additional constraints
  if (stats::coef(best_fit)["SLOPE"] <= 0.2) {
    new_start <- start
    if (cleaned_parameters$mean_resp_last > 60) {
      new_start[["IC50"]] <- min(data$log_dose, na.rm = TRUE)
    }

    new_lower <- lower
    new_lower[["SLOPE"]] <- 0.1
    new_lower[["MIN"]] <- min_lower

    new_upper <- upper
    new_upper[["SLOPE"]] <- 2.5
    new_upper[["MAX"]] <- max_upper

    best_fit <- stats::nls(
      form,
      data,
      start,
      control,
      "port",
      lower = new_lower,
      upper = new_upper
    )
  }

  best_fit
}

#' Choose fit with smallest residual standard deviation
#'
#' If there is only one non-null fit, chooses that. Otherwise, compares fits by
#' the standard deviation of their residuals and chooses the one with the
#' smallest
#'
#' This function is expected to fail in the case of two null functions
#'
#' @param fit_1,fit_2 nls fits
#' @returns nls fit
#' @noRd
choose_best_fit <- function(fit_1, fit_2) {
  # If one of them is null, no use comparing them
  if (is.null(fit_1)) {
    return(fit_2)
  }
  if (is.null(fit_2)) {
    return(fit_1)
  }

  # Otherwise, return which ever has the lowest residuals
  if (stats::sd(stats::residuals(fit_1)) < stats::sd(stats::residuals(fit_2))) {
    return(fit_1)
  }
  return(fit_2)
}
