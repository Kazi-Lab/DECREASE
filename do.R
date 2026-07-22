decrease_wrapper <- function() {
  make_fname <- function(x) paste0(x$Drug1[[1]], " & ", x$Drug2[[1]], ".RDS")

  make_final <- function() {
    # All files with & and end in .RDS. This is fragile.
    files <- fs::dir_ls(regexp = " & .*\\.RDS$")
    bind_rows(lapply(files, make_final_1))
  }

  make_final_1 <- function(path) {
    drugs <- strsplit(fs::path_ext_remove(path), " & ")[[1]]
    readRDS(path) |>
      mutate(Drug1 = drugs[1], Drug2 = drugs[2], Response = 100 - finalpred_) |>
      select(Drug1, Drug2, Conc1, Conc2, Response)
  }

  preproc <- function() {
    annot <<- read.csv(
      "ExampleData/exampleTabular/exampleTabular.csv",
      row.names = NULL
    )

    # take care of NA's and empty rows/cols
    annot <<- data.frame(lapply(annot, as.character))
    annot <<- annot[!apply(is.na(annot) | annot == "", 1, all), ] # rows with all NA
    annot <<- annot[, !apply(is.na(annot) | annot == "", 2, all)] # cols with all NA
    annot$Conc1 <<- as.numeric(annot$Conc1)
    annot$Conc2 <<- as.numeric(annot$Conc2)
    annot$Response <<- as.numeric(annot$Response)

    # convert to viability
    annot$Response <<- 100 - annot$Response
  }
  annot <<- NULL
  preproc()
  fcurve <<- FALSE
  for (k in 1:length(unique(annot$PairIndex))) {
    combiK <- annot[annot$PairIndex == unique(annot$PairIndex)[[k]], ]

    # call prediction model
    saveRDS(combiK, "annot.RDS")
    source("DeCREaSE.R")
    saveRDS(matr_Out, make_fname(combiK))
  }

  finalpred_ <<- readRDS(make_fname(data_cell_pair))
  write.csv(make_final(), "result.csv", row.names = FALSE)
}
