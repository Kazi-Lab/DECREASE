# The original package didn't really have 'function': it was more like a
# syncitium of code in a ooey gooey sac of global space. So while the old
# fit_cnmf didn't take influentPoint explicitly, it still used it.
#
# influentPoint is a boolean used to mark if an outlier was present in the data,
# which is found (in a roundabout way) using outlier_remove
#
# This argument has been renamed to be 'has outliers'
fit_cnmf <- function(data, has_outliers) {
  trainInd <- which(!is.na(data$Response))
  testInd <- which(is.na(data$Response))
  data_cell_Training <- data[trainInd, ]
  data_cell_Test <- data[testInd, ]
  cNMFpred <- do.call(
    "cbind",
    parallel::mclapply(
      1:120,
      function(i) {
        MatrTr <- reshape2::acast(
          rbind(data_cell_Training, data_cell_Test),
          Conc1 ~ Conc2,
          value.var = "Response"
        )
        MatrTr <- MatrTr +
          matrix(
            stats::runif(1, -0.001, 0.001),
            nrow = nrow(MatrTr),
            ncol = ncol(MatrTr)
          )

        if (has_outliers) {
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
  sapply(1:nrow(cNMFpred), function(i) modeest::venter(cNMFpred[i, ]))
}
