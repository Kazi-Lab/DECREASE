# One of the main problems with this function is that it (and its subfunctions
# which are sourced in) update the GLOBAL environment, making everything very
# difficult to keep track of.
#
# Another orthogonal difficulty is that the HTML elements are named non-descript
# things like 'toHide2'

# Can probably pull out a make_plot function

make_plot <- function(data, value_var, drug_1_name, drug_2_name) {
  PlotRespMatr(
    reshape2::acast(data, Conc1 ~ Conc2, value.var = value_var),
    name = "%inhibition",
    d1N = drug_1_name,
    d2N = drug_2_name
  )$pl
}
shinyServer(function(input, output, session) {
  shinyjs::runjs('$(".toHide2").hide();')
  shinyjs::disable("save_results")
  annot <<- NULL

  ########################################    On getstarted click    ########################################
  observeEvent(input$getstarted, {
    showModal(modalDialog(
      title = "DECREASE",
      fluidRow(
        column(4, fileInput("fileinp", "Upload data")),
        column(
          3,
          selectInput(
            "fileType",
            "Input file format",
            choices = c("Tabular", "Matrix")
          )
        ),
        column(
          3,
          selectInput(
            "phenotypicResponse",
            "Choose readout",
            choices = c("% cell inhibition", "% cell viability")
          )
        ),
        column(
          2,
          checkboxInput(
            "fitcurve",
            "Replace single-agent responses with fitted values"
          )
        )
      ),
      div(downloadButton(outputId = "loadExData", label = "Example data")),
      footer = tagList(modalButton("Cancel"), actionButton("start", "Start")),
    ))
    shinyjs::disable("start")
  })
  # fit curve
  observeEvent(input$fitcurve, fcurve <<- input$fitcurve)

  ###########################################    On start click    ###########################################
  observeEvent(input$start, {
    source("additionalFunctions.R")

    # change view
    shinyjs::runjs(
      '$(".hidden").removeClass("hidden");$(".toHide").hide();$(".toHide2").show();'
    )
    removeModal()

    updateSelectInput(
      session,
      "Graphs",
      choices = paste(annot$Drug1, "&", annot$Drug2)
    )
  })

  #######################################    On combination select    ########################################
  observeEvent(input$Graphs, {
    if (input$Graphs != "") {
      # selected drugs
      drugs <- unlist(strsplit(input$Graphs, " & "))
      data_cell_pair <<- annot[
        annot$Drug1 == drugs[[1]] & annot$Drug2 == drugs[[2]],
      ]
      drug_1_name <- data_cell_pair$Drug1[[1]]
      drug_2_name <- data_cell_pair$Drug2[[1]]
      pair_name <- paste0(drug_1_name, " & ", drug_2_name, ".RDS")

      if (pair_name %in% list.files()) {
        shinyjs::runjs('$(".toHide4").hide();$(".toHide3").show(); ')
        output$plotFull <- renderPlot({
          finalpred <- readRDS(pair_name)
          finalpred$finalpred_ <- 100 - finalpred$finalpred_
          make_plot(finalpred, "finalpred_", drug_1_name, drug_2_name)
        })
      } else {
        shinyjs::runjs('$(".toHide3").hide(); $(".toHide4").show();')
      }
      # plot sparce matrix to screen
      output$plotSparse <- renderPlot({
        data_cell <- data_cell_pair
        data_cell$Response <- 100 - data_cell$Response
        make_plot(data_cell, "Response", drug_1_name, drug_2_name)
      })
    }
  })

  #######################################    On start analysis click    ########################################
  observeEvent(input$startanalysis, {
    shinyjs::runjs('$(".toHide3").hide(); $(".toHide4").hide();')

    withProgress(message = 'Fitting cNMF and XGBoost', value = 0, {
      # Increment the progress bar, and update the detail text.
      incProgress(2 / 3, detail = "Running...")

      # call prediction model
      saveRDS(data_cell_pair, "annot.RDS")
      source("DeCREaSE.R")
      saveRDS(
        matr_Out,
        paste0(
          data_cell_pair$Drug1[[1]],
          " & ",
          data_cell_pair$Drug2[[1]],
          ".RDS"
        )
      )
      finalpred_ <<- matr_Out
    })

    shinyjs::runjs('$(".toHide3").show();')

    output$plotFull <- renderPlot({
      finalpred <- finalpred_
      finalpred$finalpred_ <- 100 - finalpred$finalpred_
      drug_1_name <- data_cell_pair$Drug1[[1]]
      drug_2_name <- data_cell_pair$Drug2[[1]]
      make_plot(finalpred, "finalpred_", drug_1_name, drug_2_name)
    })
    shinyjs::enable("save_results")
  })

  observeEvent(input$startanalysisall, {
    shinyjs::runjs('$(".toHide3").hide(); $(".toHide4").hide();')

    withProgress(message = 'Fitting cNMF and XGBoost', value = 0, {
      for (k in 1:length(unique(annot$PairIndex))) {
        combiK <- annot[annot$PairIndex == unique(annot$PairIndex)[[k]], ]

        # Increment the progress bar, and update the detail text.
        incProgress(
          1 / length(unique(annot$PairIndex)),
          detail = paste0(
            "Predicting ",
            combiK$Drug1[[1]],
            " & ",
            combiK$Drug2[[1]]
          )
        )

        # call prediction model
        saveRDS(combiK, "annot.RDS")
        source("DeCREaSE.R")
        saveRDS(
          matr_Out,
          paste0(combiK$Drug1[[1]], " & ", combiK$Drug2[[1]], ".RDS")
        )
      }
    })

    finalpred_ <<- readRDS(paste0(
      data_cell_pair$Drug1[[1]],
      " & ",
      data_cell_pair$Drug2[[1]],
      ".RDS"
    ))

    shinyjs::runjs('$(".toHide3").show();')

    output$plotFull <- renderPlot({
      finalpred <- finalpred_
      finalpred$finalpred_ <- 100 - finalpred$finalpred_
      drug_1_name <- data_cell_pair$Drug1[[1]]
      drug_2_name <- data_cell_pair$Drug2[[1]]
      make_plot(finalpred, "finalpred_", drug_1_name, drug_2_name)
    })
    enable("save_results")
  })

  ############################################    On file load    ############################################
  observeEvent(input$fileinp, {
    # Check that data object exists and is data frame.
    if (!is.null(input$fileinp$name)) {
      if (tools::file_ext(input$fileinp$name) %in% c("xlsx", "txt", "csv")) {
        # source file load functions
        source("getData.R")
        saveRDS(tools::file_ext(input$fileinp$name), "input$fileType")

        if (
          input$fileType == "Tabular" &&
            tools::file_ext(input$fileinp$name) %in% c("txt", "csv", "xlsx")
        ) {
          if (tools::file_ext(input$fileinp$name) == 'xlsx') {
            annot <<- openxlsx::read.xlsx(input$fileinp$datapath)
          } else if (tools::file_ext(input$fileinp$name) %in% c("txt", "csv")) {
            annot <<- read.table(
              file = input$fileinp$datapath,
              header = T,
              sep = ",",
              row.names = NULL,
              fill = T
            )
          }

          # take care of NA's and empty rows/cols
          annot <<- data.frame(lapply(annot, as.character))
          annot <<- annot[!apply(is.na(annot) | annot == "", 1, all), ] # rows with all NA
          annot <<- annot[, !apply(is.na(annot) | annot == "", 2, all)] # cols with all NA
          annot$Conc1 = as.numeric(as.character(annot$Conc1))
          annot$Conc2 = as.numeric(as.character(annot$Conc2))
          annot$Response = as.numeric(as.character(annot$Response))
        } else if (
          input$fileType == "Matrix" &&
            tools::file_ext(input$fileinp$name) %in% c("txt", "csv", "xlsx")
        ) {
          if (tools::file_ext(input$fileinp$name) == 'xlsx') {
            annot <<- openxlsx::read.xlsx(
              input$fileinp$datapath,
              colNames = F
            )
          } else if (tools::file_ext(input$fileinp$name) %in% c("txt", "csv")) {
            annot <<- read.table(
              file = input$fileinp$datapath,
              header = F,
              sep = ",",
              row.names = NULL,
              fill = T
            )
          }

          # take care of NA's and empty rows/cols
          annot <<- data.frame(lapply(annot, as.character))
          annot <<- annot[!apply(is.na(annot) | annot == "", 1, all), ] # rows with all NA
          annot <<- annot[, !apply(is.na(annot) | annot == "", 2, all)] # cols with all NA

          D1 = sum(grepl("Drug1", annot[, 1]))
          d1_ <- grep("Drug1", annot[, 1])
          annot <<- do.call(
            "rbind",
            lapply(1:D1, function(i) {
              if (i == D1) {
                annotComb <- annot[d1_[i]:nrow(annot), ]
              } else {
                annotComb <- annot[d1_[i]:(d1_[i + 1] - 1), ]
              }
              Matr = annotComb[4:nrow(annotComb), ]
              colnames(Matr) = Matr[1, ]
              Matr = Matr[-1, ]
              rownames(Matr) = Matr[, 1]
              Matr = Matr[, -1]
              Reshaped = na.omit(reshape2::melt(data.matrix(Matr)))
              colnames(Reshaped) = c("Conc1", "Conc2", "Response")
              Reshaped$ConcUnit = annotComb[3, 2]
              Reshaped$Drug1 = annotComb[1, 2]
              Reshaped$Drug2 = annotComb[2, 2]
              Reshaped$PairIndex = i
              Reshaped
            })
          )
        }
        enable("start")
        annot$Conc1 <<- as.numeric(annot$Conc1)
        annot$Conc2 <<- as.numeric(annot$Conc2)
        annot$Response <<- as.numeric(annot$Response)

        # convert to viability
        if (input$phenotypicResponse == "% cell inhibition") {
          annot$Response <<- 100 - annot$Response
        }

        # check whether enable export button
        if (
          any(sapply(1:length(unique(annot$PairIndex)), function(k) {
            combiK <- annot[
              annot$PairIndex == unique(annot$PairIndex)[[k]],
            ]
            file.exists(paste0(
              combiK$Drug1[[1]],
              " & ",
              combiK$Drug2[[1]],
              ".RDS"
            ))
          }))
        ) {
          enable("save_results")
        }
      }
    }
  })

  # Upload example data
  output$loadExData <- downloadHandler(
    "ExampleData.zip",
    \(file) file.copy("ExampleData.zip", file)
  )

  observeEvent(input$save_results, {
    pairsCalculated_ <- as.character(na.omit(sapply(
      1:length(unique(annot$PairIndex)),
      function(k) {
        combiK <- annot[annot$PairIndex == unique(annot$PairIndex)[[k]], ]
        if (
          file.exists(
            paste0(combiK$Drug1[[1]], " & ", combiK$Drug2[[1]], ".RDS")
          )
        ) {
          paste(combiK$Drug1[[1]], "&", combiK$Drug2[[1]])
        } else {
          NA
        }
      }
    )))

    showModal(modalDialog(
      title = "Choose options",

      fluidRow(
        column(
          width = 8,
          offset = 1,
          selectizeInput(
            "toExportDrugs",
            "Choose drug pairs",
            choices = pairsCalculated_,
            multiple = T,
            selected = pairsCalculated_[1:length(pairsCalculated_)],
            width = "100%"
          )
        ),
        column(
          width = 2,
          downloadButton("downloadPDF", label = "Download *.pdf")
        )
      ),

      fluidRow(
        column(
          2,
          offset = 1,
          downloadButton("downloadMatrices", label = "Download (.xlsx)")
        ),
        column(
          2,
          downloadButton("downloadMatrices2", label = "Download (.csv)")
        ),
        column(
          3,
          offset = 1,
          selectInput(
            "compatib",
            label = "Compatible with:",
            choices = c("SynergyFinder", "Combenefit"),
            selected = "SynegyFinder"
          )
        )
      ),
    ))

    output$downloadPDF <- downloadHandler(
      \() paste0("result_", Sys.Date(), ".pdf"),
      \(file) {
        pdf("fileOut.pdf", width = 12)
        for (k in 1:length(input$toExportDrugs)) {
          drugs <- unlist(strsplit(input$toExportDrugs[k], " & "))
          combiK <- annot[annot$Drug1 == drugs[1] & annot$Drug2 == drugs[2], ]
          combiApprox <- readRDS(
            paste0(combiK$Drug1[[1]], " & ", combiK$Drug2[[1]], ".RDS")
          )

          drug_1_name <- data_cell_pair$Drug1[[1]]
          drug_2_name <- data_cell_pair$Drug2[[1]]
          combiApprox$Response <- 100 - combiApprox$Response
          combiApprox$finalpred_ <- 100 - combiApprox$finalpred_
          gridExtra::grid.arrange(
            make_plot(combiApprox, "Response", drug_1_name, drug_2_name),
            make_plot(combiApprox, "finalpred_", drug_1_name, drug_2_name),
            ncol = 2
          )
        }
        dev.off()

        file.copy("fileOut.pdf", file)
      },
    )

    # create output compatible with synergyfinder
    synergyMakeFinal <- function() {
      do.call(
        "rbind",
        lapply(1:length(input$toExportDrugs), function(k) {
          combiK <- annot[
            annot$Drug1 == unlist(strsplit(input$toExportDrugs[k], " & "))[1] &
              annot$Drug2 == unlist(strsplit(input$toExportDrugs[k], " & "))[2],
          ]
          combiApprox <- readRDS(paste0(
            combiK$Drug1[[1]],
            " & ",
            combiK$Drug2[[1]],
            ".RDS"
          ))
          combiApprox$Drug1 = combiK$Drug1[[1]]
          combiApprox$Drug2 = combiK$Drug2[[1]]
          combiApprox$PairIndex = k
          combiApprox$Response = combiApprox$finalpred_
          combiApprox$ConcUnit = combiK$ConcUnit[[1]]
          combiApprox$Response = 100 - combiApprox$Response
          combiApprox
        })
      )
    }

    out_cols <- c(
      "PairIndex",
      "Drug1",
      "Drug2",
      "Conc1",
      "Conc2",
      "Response",
      "ConcUnit"
    )
    output$downloadMatrices <- downloadHandler(
      \() paste0("result_", Sys.Date(), ".xlsx"),
      \(file) {
        out <- synergyMakeFinal()
        openxlsx::write.xlsx(out[, out_cols], "fileOut.xlsx")
        file.copy("fileOut.xlsx", file)
      },
    )
    output$downloadMatrices2 <- downloadHandler(
      \() paste0("result_", Sys.Date(), ".csv"),
      \(file) {
        out <- synergyMakeFinal()
        write.csv(out[, out_cols], "fileOut.csv", row.names = FALSE)
        file.copy("fileOut.csv", file)
      }
    )
  })

  session$onSessionEnded(\() {
    stopApp()
    quit("no")
  })
})
