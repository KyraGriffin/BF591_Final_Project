## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)


deseq_choices <-
    c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")


# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = bs_theme(bootswatch = "minty"),
    # Application title
    mainPanel(tabsetPanel(
      id = "tabsetPanelID",
      type = "tabs",
      tabPanel("Tab1", tabsetPanel(
        tabPanel("SubPanelA1"), tabPanel("SubPanelA2")
      )),
      tabPanel("Tab2", tabsetPanel(
        tabPanel("SubPanelB1"), tabPanel("SubPanelB2")
      )),
      tabPanel("Tab3", tabsetPanel(
        tabPanel("SubPanelC1"), tabPanel("SubPanelC2")
      ))
    ))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
}

# Run the application
shinyApp(ui = ui, server = server)
