lapply(
  c(
    "scales",
    "drc",
    "reshape2",
    "xgboost",
    "parallel",
    "openxlsx",
    "dplyr",
    "NNLM",
    "modeest",
    "mlrMBO",
    "lhs",
    "caret"
  ),
  require,
  character.only = !0
)

outlier_remove <- function(xTmp, iqr_ = 1.5) {
  qq <- unname(quantile(xTmp, probs = c(.25, .75), na.rm = T))
  outlier_detector <- iqr_ * IQR(xTmp, na.rm = T)
  xTmp < (qq[1] - outlier_detector) | xTmp > (qq[2] + outlier_detector)
}

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

#' @param df data.frame with Conc1, Conc2, and Response
convert_to_matrix <- function(df) {
  reshape2::acast(df, Conc1 ~ Conc2, value.var = "Response")
}

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

gdes <- function() {
  generateDesign(
    n = 120,
    par.set = makeParamSet(
      makeNumericVectorParam("logNtree", len = 1, lower = 4, upper = 9),
      makeNumericVectorParam("lambda", len = 1, lower = 0, upper = 3),
      makeNumericVectorParam("alpha", len = 1, lower = 0, upper = 3),
      makeIntegerVectorParam("maxdepth", len = 1, lower = 1, upper = 6),
      makeNumericVectorParam("subsample", len = 1, lower = .4, upper = 1),
      makeNumericVectorParam(
        "colsample_bytree",
        len = 1,
        lower = .4,
        upper = 1
      ),
      makeNumericVectorParam("eta", len = 1, lower = .001, upper = .1)
    ),
    fun = lhs::randomLHS
  )
}
#' @param path character path to .csv file containing columns Conc1, Conc2, and
#'   Response
#' @param is_viability boolean. Is the Response column viability? If FALSE,
#'   converted to viability.
#' @param use_fitted_single_agent_values boolean. Should the fitted single-agent values be
#'   used? If FALSE, use raw values.
dc <- function(path, is_viability, use_fitted_single_agent_values) {
  data_cell <- read_data(path) |>
    preprocess_data(is_viability)
  set.seed(42)
  influentPoint <- NULL # for now
  MatrTr <- convert_to_matrix(data_cell)
  MatrTr <- make_0_drug_highest(MatrTr)

  # LOOK OUT! MatrTr is now inhibition
  MatrTr <- 100 - MatrTr

  bliss.mat <- calculate_bliss(MatrTr)

  # fit single-agent
  # TODO: Compat function that takes D1Len (and D2Len) and converts it to
  # dataframe with dose and response cols
  #
  # Ensure response is correct! Not sure if should be inhibition or viability
  D1Len <- MatrTr[, 1]
  names(D1Len)[1] <- 1e-6
  D2Len <- MatrTr[1, ]
  names(D2Len)[1] <- 1e-6
  d1 <- tryCatch(predict(fit_dose_response(D1Len)), error = \(e) D1Len)
  d2 <- tryCatch(predict(fit_dose_response(D2Len)), error = \(e) D2Len)

  # new fix half curve
  MatrTrCopy <- MatrTr
  if (use_fitted_single_agent_values) {
    MatrTr[, 1] <- d1
    MatrTr[1, ] <- d2
  } else {
    MatrTr[, 1] <- (MatrTr[, 1] + d1) / 2
    MatrTr[1, ] <- (MatrTr[1, ] + d2) / 2
  }

  # single-agent deviations
  devD1 <- abs(d1 - (MatrTr[, 1]))
  devD2 <- abs(d2 - (MatrTr[1, ]))
  dev_ <- abs(bliss.mat - MatrTr)
  MatrOutl <- matrix(!1, nrow = nrow(MatrTr), ncol = ncol(MatrTr))

  # check deviations with Bliss
  MatrOutl[-1, -1] <- outlier_remove(
    abs(dev_[-1, -1] - median(dev_[-1, -1], na.rm = T)),
    iqr_ = 5
  ) &
    (dev_[-1, -1] > 25)
  MatrOutl[, 1] <- as.logical(colSums(MatrOutl, na.rm = T)) &
    devD1 > 10 |
    devD1 > 15
  MatrOutl[1, ] <- as.logical(rowSums(MatrOutl, na.rm = T)) &
    devD2 > 10 |
    devD1 > 15

  # remove possible outliers
  MatrOutl[is.na(MatrOutl)] <- !1
  MatrTr[MatrOutl] <- NA
  if (any(MatrOutl)) {
    influentPoint <- T
    warning("Possible outliers were removed")
  }

  matr_Out <- reshape2::melt(100 - MatrTr)
  colnames(matr_Out) <- c("Conc1", "Conc2", "Response")
  matr_Out <- dplyr::arrange(matr_Out, Conc1, Conc2)
  matr_OutCopy <- reshape2::melt(100 - MatrTrCopy)
  colnames(matr_OutCopy) <- c("Conc1", "Conc2", "Response")
  matr_OutCopy <- dplyr::arrange(matr_OutCopy, Conc1, Conc2)

  # fill single-agent response columns (check the design section of manuscript)
  matr_Out$R1 <- sapply(1:nrow(matr_Out), function(i) {
    matr_Out$Response[
      matr_Out$Conc1 == matr_Out[i, ]$Conc1 & matr_Out$Conc2 == 0
    ]
  })
  matr_Out$R2 <- sapply(1:nrow(matr_Out), function(i) {
    matr_Out$Response[
      matr_Out$Conc2 == matr_Out[i, ]$Conc2 & matr_Out$Conc1 == 0
    ]
  })

  # training and test
  trainInd <- which(!is.na(matr_Out$Response))
  testInd <- which(is.na(matr_Out$Response))
  data_cell_Training <- matr_Out[trainInd, ]
  data_cell_Test <- matr_Out[testInd, ]

  ######################################################################################################################
  ################################################     fit cNMF      ###################################################

  cNMFpred <- do.call(
    "cbind",
    mclapply(
      1:120,
      function(i) {
        MatrTr <- reshape2::acast(
          rbind(data_cell_Training, data_cell_Test),
          Conc1 ~ Conc2,
          value.var = "Response"
        )
        MatrTr <- MatrTr +
          matrix(
            runif(1, -0.001, 0.001),
            nrow = nrow(MatrTr),
            ncol = ncol(MatrTr)
          )

        if (length(influentPoint) != 0) {
          if (is.na(MatrTr[1, 1])) {
            MatrTr[1, 1] <- 100
          }
          if (is.na(MatrTr[nrow(MatrTr), 1])) {
            MatrTr[nrow(MatrTr), 1] <- MatrTr[nrow(MatrTr) - 1, 1]
          }
          if (is.na(MatrTr[1, ncol(MatrTr)])) {
            MatrTr[1, ncol(MatrTr)] <- MatrTr[1, ncol(MatrTr) - 1]
          }
          MatrTr[, 1] <- zoo::na.approx(MatrTr[, 1], rule = 2)
          MatrTr[1, ] <- zoo::na.approx(MatrTr[1, ], rule = 2)
        }

        nsclc2.nmf <- NNLM::nnmf(
          as.matrix(MatrTr),
          sample(2:3, 1),
          verbose = F,
          beta = sample(seq(.1, 1, .01), 3),
          alpha = sample(seq(.1, 1, .01), 3),
          max.iter = 500L,
          loss = "mse",
          check.k = F
        )
        nsclc2.hat.nmf <- with(nsclc2.nmf, W %*% H)

        matr_Out <- reshape2::melt(nsclc2.hat.nmf)
        colnames(matr_Out) <- c("Conc1", "Conc2", "Response")
        dplyr::arrange(matr_Out, Conc1, Conc2)$Response
      },
      mc.cores = 8
    )
  )

  # prepare predictions
  if (sum(colSums(cNMFpred == 0) == 0) > 1) {
    cNMFpred <- cNMFpred[, colSums(cNMFpred == 0) == 0]
  }

  cNMFpred <- sapply(1:nrow(cNMFpred), function(i) {
    modeest::venter(cNMFpred[i, ])
  })

  ### parameter set

  obj.fun.err <- function(x) {
    logNtree <- x[1][[1]]
    lambda <- x[2][[1]]
    alpha <- x[3][[1]]
    maxdepth <- x[4][[1]]
    subsample <- x[5][[1]]
    colsample_bytree <- x[6][[1]]
    eta <- x[7][[1]]
    MAE_ <- 0

    # repeated CV
    for (repCv in 1:2) {
      flds <- caret::createFolds(
        data_cell_Training$Response,
        k = 3,
        list = T,
        returnTrain = F
      )
      for (k in 1:length(flds)) {
        testData <- data_cell_Training[flds[[k]], ]
        trainData <- data_cell_Training[-flds[[k]], ]
        fit <- xgboost(
          data = data.matrix(
            trainData[,
              c(
                "R1",
                "R2",
                "Conc1",
                "Conc2"
              )
            ]
          ),
          label = trainData$Response,
          verbose = F,
          nrounds = round(2**logNtree),
          nthread = 8,
          params = list(
            objective = "reg:linear",
            max.depth = maxdepth,
            eta = eta,
            lambda = lambda,
            alpha = alpha,
            subsample = subsample,
            colsample_bytree = colsample_bytree
          )
        )

        ypred <- predict(
          fit,
          as.matrix(testData[, c("R1", "R2", "Conc1", "Conc2")])
        )
        MAE_ <- MAE_ + mean(abs(ypred - testData$Response), na.rm = T)
      }
    }
    MAE_
  }
  des <- gdes()
  gc(T)
  des$y <- apply(des, 1, obj.fun.err)

  models <- des[order(des$y), ]
  # Fit with 5 models with best parameters
  XGBoostpred <- do.call(
    "cbind",
    lapply(1:5, function(i) {
      fit <- xgboost(
        as.matrix(data_cell_Training[, c("R1", "R2", "Conc1", "Conc2")]),
        label = data_cell_Training$Response,
        verbose = F,
        nrounds = round(2**models[i, "logNtree"]),
        nthread = 8,
        params = list(
          objective = "reg:linear",
          max.depth = models[i, "maxdepth"],
          eta = models[i, "eta"],
          lambda = models[i, "lambda"],
          alpha = models[i, "alpha"],
          subsample = models[i, "subsample"],
          colsample_bytree = models[i, "colsample_bytree"]
        )
      )
      predict(fit, as.matrix(matr_Out[, c("R1", "R2", "Conc1", "Conc2")]))
    })
  )
  XGBoostpred <- sapply(1:nrow(XGBoostpred), function(i) {
    modeest::venter(XGBoostpred[i, ])
  })

  # final prediction
  finalpred_ <- (XGBoostpred + cNMFpred) / 2

  matr_Out$cNMFpred <- cNMFpred
  matr_Out$XGBoostpred <- XGBoostpred
  matr_Out$finalpred_ <- finalpred_
  matr_Out$finalpred_[trainInd] <- matr_Out$Response[trainInd]
  if (!use_fitted_single_agent_values) {
    matr_Out[
      matr_Out$Conc1 == 0 | matr_Out$Conc2 == 0,
      "finalpred_"
    ] <- matr_OutCopy[
      matr_OutCopy$Conc1 == 0 | matr_OutCopy$Conc2 == 0,
      "Response"
    ]
  }
}
