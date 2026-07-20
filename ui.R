library(shiny)
library(shinyjs)
library(gridExtra)

shinyUI(navbarPage(
  useShinyjs(),
  tags$head(tags$script(src = "plotsurface.js")),

  div(
    class = "toHide2",
    column(
      2,
      offset = 0,
      div(
        class = "hidden",
        fluidRow(
          div(
            class = "selectWrap",
            selectInput("Graphs", "Select combination:", choices = "")
          )
        ),
        fluidRow(actionButton("save_results", "Export results"))
      )
    ),
    column(4, offset = 1, plotOutput('plotSparse')),
    column(
      4,
      offset = 1,
      div(class = "toHide3", plotOutput('plotFull')),
      div(
        class = "hidden toHide4",
        div(
          fluidRow(HTML('Predict current')),
          fluidRow(actionButton("startanalysis", "Start"))
        ),
        div(
          fluidRow(HTML('Predict all')),
          fluidRow(actionButton("startanalysisall", "Start all"))
        )
      )
    ),
  ),
  fluidRow(
    div(
      class = "toHide",
      column(2, offset = 4, actionButton("getstarted", "Get started"))
    )
  ),
))
