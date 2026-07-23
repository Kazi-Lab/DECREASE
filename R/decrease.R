#' Impute missing values from dose-response experiments
#'
#' @param data data.frame containing columns Conc1, Conc2, and Response
#' @param is_viability boolean. Is the Response column viability? If FALSE,
#'   converted to viability.
#' @param use_fitted_single_agent_values boolean. Should the fitted single-agent
#'   values be used? If FALSE, use raw values.
#' @export
decrease <- function(data, is_viability, use_fitted_single_agent_values) {
  # TODO: implement is_viability
  set.seed(42)

  data <- preprocess_data(data, is_viability)

  mat <- convert_to_matrix(data)
  mat <- make_0_drug_highest(mat)

  # LOOK OUT! mat is now inhibition
  mat <- 100 - mat

  bliss.mat <- calculate_bliss(mat)

  get_fitted_responses <- function(responses) {
    names(responses)[1] <- 1e-6
    responses <- dplyr::as_tibble(responses, rownames = "dose")
    colnames(responses) <- c("dose", "response")
    tryCatch(predict(fit_dose_responses(responses)), \(x) responses)
  }

  outliers <- detect_outliers(
    mat,
    use_fitted_single_agent_values,
    get_fitted_responses(mat[, 1]),
    get_fitted_responses(mat[1, ])
  )

  if (any(outliers, na.rm = TRUE)) {
    warning("Possible outliers were removed")
  }

  # LOOK OUT! out is viability now
  out <- reshape2::melt(100 - mat)
  colnames(out) <- c("Conc1", "Conc2", "Response")
  out <- dplyr::arrange(out, Conc1, Conc2)
  og <- out

  # fill single-agent response columns (check the design section of manuscript)
  out$R1 <- sapply(
    1:nrow(out),
    \(i) out$Response[out$Conc1 == out[i, ]$Conc1 & out$Conc2 == 0]
  )
  out$R2 <- sapply(
    1:nrow(out),
    \(i) out$Response[out$Conc2 == out[i, ]$Conc2 & out$Conc1 == 0]
  )

  trainInd <- which(!is.na(out$Response))
  out$cNMFpred <- fit_cnmf(out, influentPoint)
  out$XGBoostpred <- fit_xgboost(out)
  out$finalpred_ <- (out$XGBoostpred + out$cNMFpred) / 2
  out$finalpred_[trainInd] <- out$Response[trainInd]

  if (!use_fitted_single_agent_values) {
    out[out$Conc1 == 0 | out$Conc2 == 0, "finalpred_"] <-
      og[og$Conc1 == 0 | og$Conc2 == 0, "Response"]
  }

  out
}
