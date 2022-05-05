## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(vtable)
library(tidyverse)
library(bslib)
library(ggplot2)
library(colourpicker)

dataset_choice <- c("Neurologically normal", "Huntington's Disease")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Final Project"),
  tabsetPanel(
    tabPanel(
      "Samples",
      sidebarLayout(
        sidebarPanel(fileInput(
          "sample_metadata",
          label = "Load sample information matrix",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          ),
          placeholder = "sample_metadata.csv"
        ), radioButtons(
          inputId = "diagnosis_samp",
          label = "Choose the Diagnosis",
          choices = dataset_choice,
          selected = "Neurologically normal"
        ),
        actionButton(inputId = "button", label = "Submit")),
        mainPanel(
          tabsetPanel(
            id = "tabsetPanelID",
            type = "tabs",
            tabsetPanel(
              tabPanel("Summary", tableOutput("samp_sum_table")),
              tabPanel("Table"),
              tabPanel("Plots")
            )
          ),
        )
      )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 35 * 1024^2)

  #' load_Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is
  #' still a "function", but it will take no arguments. The `reactive({})` bit
  #' says "if any of my inputs (as in, input$...) are changed, run me again".
  #' This is useful when a user clicks a new button or loads a new file. In
  #' our case, look for the uploaded file's datapath argument and load it with
  #' read.csv. Return this data frame in the normal return() style.
  load_sample_data <- reactive({
    inFile <- input$sample_metadata
    
    if (is.null(inFile)) {
      return(NULL)
    }
    data <- read.csv(file = inFile$datapath, header=TRUE, stringsAsFactors = TRUE)
    
    return(data)
  })

  #' Draw summary table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param diagnosis Diagnosis of sample row from the radio button input.
  #'
  #' @return Data frame
  #'
  #' @details I would suggest the function
  #' `formatC()`
  #'
  #' @examples draw_sum_table(sample_metadata, "Neurologically normal")
  draw_sum_table <- function(dataf, diagnosis) {
    df_out <- dataf %>% 
      dplyr::filter(Diagnosis == diagnosis)
    
    if(diagnosis == 'Neurologically normal'){
      df_out <- dataf %>% 
        dplyr::select(-c(age_of_onset, cag, Duration, h.v_cortical_score,h.v_striatal_score,vonsattel_grade))
      
      column_names <- colnames(df_out)
      type_col <- sapply(df_out,class)
      distinct_val <- sapply(df_out, function(x) toString(unique(x)[1:5]))
      mean_sd_col <- sapply(df_out,is.factor)
      
    } else{
      
      column_names <- colnames(df_out)
      type_col <- sapply(df_out,class)
      distinct_val <- sapply(df_out, function(x) toString(unique(x)))
      mean_sd_col <- sapply(df_out,is.factor)
    }

    df <- data.frame(column_names, type_col, distinct_val, mean_sd_col)
    
    #colnames(df) = c("Column Name",	"Type", "Mean (sd) or Distinct Values", "temp vals")
    return(df)
  }

  # Sample Summary Table
  output$samp_sum_table <- renderTable(
    {
      req(input$sample_metadata)
      table <- load_sample_data()
      return(draw_sum_table(dataf = table, diagnosis = input$diagnosis_samp))
    },
    striped = T
  )
}

# Run the application
shinyApp(ui = ui, server = server)
