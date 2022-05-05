## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(DT)
library(tidyverse)
library(bslib)
library(ggplot2)
library(colourpicker)

dataset_choice <- c("Neurologically normal", "Huntington's Disease")
sample_y_choice <- c("age_of_death", "AvgSpotLen","Bases","Bytes","mrna.seq_reads","pmi","rin",
                     "age_of_onset","cag","Duration","h.v_cortical_score","h.v_striatal_score","vonsattel_grade")

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
              tabPanel("Table", DT::dataTableOutput("sample_DT")),
              tabPanel("Plots", radioButtons(
                inputId = "samp_y_axis",
                inline = TRUE,
                label = "Choose the y-axis",
                choices = sample_y_choice,
                selected = "age_of_death"
              ),plotOutput("sample_plot"))
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
  
  #' Draw sample plot
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
  draw_sample_plot <- function(dataf, y_val) {
    df_out <- dataf 
    hd <- c("age_of_onset","cag","Duration","h.v_cortical_score","h.v_striatal_score","vonsattel_grade")
    
    plot_data <- df_out %>% dplyr::select(where(is.numeric))
    plot_data <- plot_data %>%
      add_column(Diagnosis = df_out$Diagnosis)
    
    if(y_val %in% hd){
      plot_data <- plot_data %>% dplyr::filter(Diagnosis == "Huntington's Disease")
    }
    
    p <- ggplot(plot_data, aes(x = Diagnosis, y = !!sym(y_val), fill=Diagnosis))+ # aes(x = !!sym(x_name), y = !!sym(y_name)
      geom_violin()
      theme_bw() + 
      theme(legend.position = "bottom")
    
    return(p)
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
  
  output$sample_DT <- DT::renderDataTable({
    req(input$sample_metadata)
    table <- load_sample_data()
    table <- table %>% 
      dplyr::filter(Diagnosis == input$diagnosis_samp)
    if(input$diagnosis_samp == 'Neurologically normal'){
      table <- table %>%
        dplyr::select(-c(age_of_onset, cag, Duration, h.v_cortical_score,h.v_striatal_score,vonsattel_grade)) 
    }
    DT::datatable(table, options = list(orderClasses = TRUE))
  })
  
  output$sample_plot <- renderPlot(
    {
      req(input$sample_metadata)
      table <- load_sample_data()
      if(!is.null(table)){
        draw_sample_plot(table, 
                   input$samp_y_axis)}})
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
