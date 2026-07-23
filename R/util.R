read_data <- function(path) read.csv(path, row.names = NULL)

#' Check data integrity and format for analysis
#'
#' This function checks for the presence of columns as well as if they are of
#' appropriat type, and converts Response to viability if is_viability is FALSE.
#'
#' @param data data.frame containing columns Conc1, Conc2, and Response, all of
#'   which should be numeric
#' @param is_viability boolean. Is the Response column measuring viability? If
#'   not, assumes measuring inhibition, and will be converted to viability.
preprocess_data <- function(data, is_viability) {
  # Unlike in the original, this function does not check for rows and/or cols
  # that are entirely NA.
  required_names <- c("Conc1", "Conc2", "Response")
  stopifnot(
    is.data.frame(data),
    all(required_names %in% names(data)),
    is.numeric(data$Conc1),
    is.numeric(data$Conc2),
    is.numeric(data$Response)
  )

  if (!is_viability) {
    message("Converting Response to viability")
    data$Response <- 100 - data$Response
  }
  data
}

#' Ensure 0 drug condition has highest viability
#'
#' If the 0 drug condition is lower than either of its lowest-concentration
#' single agent drugs, set its viability to 100
#'
#' @param mat a matrix of viabilities, where the top left entry is the 0, 0 drug
#'   condition
make_0_drug_highest <- function(mat) {
  if (mat[1, 1] < max(mat[2, 1], mat[1, 2])) {
    mat[1, 1] <- 100
  }

  mat
}

#' Convert dose-response data.frame into a matrix
#'
#' @param df data.frame with Conc1, Conc2, and Response
convert_to_matrix <- function(df) {
  # This function exists largely so if I decide to drop reshape2 and use dplyr
  # instead, I don't have to go hunting as much.
  reshape2::acast(df, Conc1 ~ Conc2, value.var = "Response")
}
