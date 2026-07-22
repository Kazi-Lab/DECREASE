fit_xgboost <- function(data) {
  trainInd <- which(!is.na(data$Response))
  data_cell_Training <- data[trainInd, ]

  ### parameter set
  gdes <- function() {
    ParamHelpers::generateDesign(
      n = 120,
      par.set = ParamHelpers::makeParamSet(
        ParamHelpers::makeNumericVectorParam(
          "logNtree",
          len = 1,
          lower = 4,
          upper = 9
        ),
        ParamHelpers::makeNumericVectorParam(
          "lambda",
          len = 1,
          lower = 0,
          upper = 3
        ),
        ParamHelpers::makeNumericVectorParam(
          "alpha",
          len = 1,
          lower = 0,
          upper = 3
        ),
        ParamHelpers::makeIntegerVectorParam(
          "maxdepth",
          len = 1,
          lower = 1,
          upper = 6
        ),
        ParamHelpers::makeNumericVectorParam(
          "subsample",
          len = 1,
          lower = .4,
          upper = 1
        ),
        ParamHelpers::makeNumericVectorParam(
          "colsample_bytree",
          len = 1,
          lower = .4,
          upper = 1
        ),
        ParamHelpers::makeNumericVectorParam(
          "eta",
          len = 1,
          lower = .001,
          upper = .1
        )
      ),
      fun = lhs::randomLHS
    )
  }

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
        fit <- xgboost::xgboost(
          data = data.matrix(trainData[, c("R1", "R2", "Conc1", "Conc2")]),
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
      fit <- xgboost::xgboost(
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
      predict(fit, as.matrix(data[, c("R1", "R2", "Conc1", "Conc2")]))
    })
  )
  XGBoostpred <- sapply(
    1:nrow(XGBoostpred),
    function(i) modeest::venter(XGBoostpred[i, ])
  )
  XGBoostpred
}
