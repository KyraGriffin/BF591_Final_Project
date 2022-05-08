## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(DT)
library(plotly)
library(heatmaply)
library(tidyverse)
library(bslib)
library(ggplot2)
library(plotly)
library(colourpicker)
#install.packages('shinyHeatmaply')

dataset_choice <- c("Neurologically normal", "Huntington's Disease")
sample_y_choice <- c("age_of_death", "AvgSpotLen","Bases","Bytes",
                     "mrna.seq_reads","pmi","rin","age_of_onset","cag",
                     "Duration","h.v_cortical_score","h.v_striatal_score",
                     "vonsattel_grade")

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
    ),
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(fileInput(
          "count_file",
          label = "Load Normalized counts matrix",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          ),
          placeholder = "GSE64810_mlhd_DESeq2_norm_counts_adjust.csv"
        ), 
        sliderInput(
          "variance_slider",
          "Select the percentile of variance:",
          min = 0,
          max = 100,
          value = 50,
          step = 1
        ),
        sliderInput(
          "non_zero_slider",
          "Select the number of the samples that are non-zero:",
          min = 0,
          max = 100,
          value = 0,
          step = 1
        ),
        actionButton(inputId = "count_button", label = "Submit")),
        mainPanel(
          tabsetPanel(
            id = "tabsetPanelID",
            type = "tabs",
            tabsetPanel(
              tabPanel("Summary", br(), tableOutput("count_sum_table")), 
              tabPanel("Diagnostic Scatter Plot", br(),
                       splitLayout(cellWidths = c("50%", "50%"), 
                                   plotOutput("median_vs_var_plot"), 
                                   plotOutput("median_vs_zeros"))), 
              tabPanel("Clustered Heatmap",  br(),
                       plotOutput("count_heatmap")),
              tabPanel("PCA",
                       h5("Select PC vs PC or Beeswarm to plot either two principle components or 
                          all the principle components together in a Beeswarm Plot:"),
                       radioButtons("pc_button", "",
                                    c("PC vs PC", "Beeswarm Plot")),
                       h5("If you choose PC vs PC please enter the numbers of the two principle components you want to plot below:"),
                       textInput("PC1", label = h6("First PC"), value = "...", width = 50),
                       textInput("PC2", label = h6("Second PC"), value = "...", width = 70),
                       plotOutput("count_pca"))
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
  
  ############### Load Data ####################
  
  #' Load Sammple Data
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
  
  #' Load Counts Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is 
  #' still a "function", but it will take no arguments. The `reactive({})` bit 
  #' says "if any of my inputs (as in, input$...) are changed, run me again". 
  #' This is useful when a user clicks a new button or loads a new file. In 
  #' our case, look for the uploaded file's datapath argument and load it with 
  #' read.csv. Return this data frame in the normal return() style.
  load_count_data <- reactive({
    df <- read_tsv(input$count_file$datapath, col_names = TRUE, show_col_types = FALSE)
    colnames(df)[1] <- "gene"
    return(df)
  })

  ############### Functions ####################
  
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
  
  #' Draw count summary table
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
  draw_count_sum_table <- function(count_data, var_val, zero_val) {
    #number of samples
    num_samples <- NCOL(count_data)
    #total number of genes
    total_genes <- length(count_data$gene)
    
    # Filter data based on variance percentile and nonzeros
    nonzero_genes <- rowSums(count_data[-1] == 0)>= zero_val
    nonzero_counts <- count_data[nonzero_genes,]
    zero_counts <- filter(count_data, !gene %in% nonzero_counts$gene)
    nonzero_counts$variance <- apply(nonzero_counts[,-c(1)], 1, var)
    
    percent = quantile(nonzero_counts$variance, prob = var_val / 100)
    
    var_counts <- dplyr::filter(nonzero_counts, variance >= percent) %>% 
      select(-variance)
    
    #number and % of genes passing current filter
    num_pass <- length(var_counts$gene)
    n1 <- round((num_pass / total_genes) * 100, 2)
    pass <- c(num_pass, n1)
    num_pass <- toString(pass, width = 30)
    
    #number and % of genes not passing current filter
    np <- length(zero_counts$gene) + (length(nonzero_counts$gene) - length(var_counts$gene))
    n2 <- round(((np / total_genes) * 100), 2)
    not_pass <- c(np,n2)
    
    num_not_pass <- toString(not_pass, width = 30)

    df <- data.frame(num_samples, total_genes, num_pass, num_not_pass)
    colnames(df) = c("Number of Samples",	"Total Number of Genes", 
                     "Number , % of genes passing the filter", 
                     "Number , % of genes not passing the filter")
    return(df)
  }
  
  
  #' Draw count summary table
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
  filter_count_data <- function(count_data, var_val, zero_val) {

    # Filter data based on variance percentile and nonzeros
    nonzero_genes <- rowSums(count_data[-1] == 0)>= zero_val
    nonzero_counts <- count_data[nonzero_genes,]
    
    zero_counts <- filter(count_data, !gene %in% nonzero_counts$gene)
    nonzero_counts$variance <- apply(nonzero_counts[,-c(1)], 1, var)
    
    percent = quantile(nonzero_counts$variance, prob = var_val / 100)
    var_counts <- dplyr::filter(nonzero_counts, variance >= percent) %>% 
      select(-variance)
    not_var_counts <- filter(nonzero_counts, !gene %in% var_counts$gene) %>% 
      select(-variance)
    
    not_filtered <- rbind(zero_counts, not_var_counts) %>% mutate(volcano = "Not Filtered")
    filtered <- var_counts %>% mutate(volcano = "Filtered")
    
    data <- as.data.frame(rbind(filtered, not_filtered))
    

    return(data)
  }
  
  var_volcano_plot <-function(data,x_lab, y_lab, title) {
    dt <- data %>% select(-volcano)
    medians <- apply(dt[,-c(1)], 1, median)
    variances <- apply(dt[,-c(1)], 1, var)
    
    plot_data <- tibble::tibble(median=medians, variance=variances, volcano = data$volcano)
    plot_data$rank <- rank(plot_data$median)
    
    p <- ggplot(plot_data, aes(x=rank, y=variance)) +
      geom_point(aes(color = volcano)) + 
      scale_color_manual(values = c("#F8766D","#619CFF")) + 
      theme_bw() +
      ggplot2::xlab(x_lab) +
      ggplot2::ylab(y_lab) +
      ggtitle(title)
    theme(legend.position = "bottom")
    
    return(p)
  }
  
  
  zero_volcano_plot <-function(data,x_lab, y_lab, title) {
    dt <- data %>% select(-volcano)
    medians <- apply(dt[,-c(1)], 1, median)
    num_zero <- rowSums(dt[,-c(1)] == 0)
    
    plot_data <- tibble::tibble(median=medians, zeros=num_zero, volcano = data$volcano)
    plot_data$rank <- rank(plot_data$median)
    
    p <- ggplot(plot_data, aes(x=rank, y=zeros)) +
      geom_point(aes(color = volcano)) + 
      scale_color_manual(values = c("#F8766D","#619CFF")) + 
      theme_bw() +
      ggplot2::xlab(x_lab) +
      ggplot2::ylab(y_lab) +
      ggtitle(title)
      theme(legend.position = "bottom")
          
      return(p)
  }
  
  plot_heatmap <- function(data_filtered, title) {
    d_filter <- data_filtered %>% dplyr::filter(volcano != "Not Filtered") %>% 
      select(-volcano) %>% 
      column_to_rownames(var = "gene")
    
    print(select_if(d_filter, is.numeric))
    
    p <- heatmap(as.matrix(d_filter), 
                 margins = c(5,5), 
                 main = title,
                 cexRow =0.7,
                 cexCol = 0.7)
    #p <- heatmaply::heatmaply(as.matrix(d_filter))
    
    return(p)
  }

  ############### Output ####################
  
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
  
  output$count_sum_table <- renderTable(
    {
      req(input$count_file)
      table <- load_count_data()
      if(!is.null(table)){
        draw_count_sum_table(table, 
                       input$variance_slider,
                       input$non_zero_slider)}}, spacing = "m",
    bordered = T)
  
  output$median_vs_var_plot <- renderPlot(
    {
      req(input$count_file)
      table <- load_count_data()
      if(!is.null(table)){
        data <- filter_count_data(table,
                                  input$variance_slider,
                                  input$non_zero_slider)
        var_volcano_plot(data,
                         "Rank(Median)",
                         "Variance",
                         "Median Read Counts and Variability Over All Genes")}})
  
  
  output$median_vs_zeros <- renderPlot(
    {
      req(input$count_file)
      table <- load_count_data()
      if(!is.null(table)){
      data <- filter_count_data(table,
                                input$variance_slider,
                                input$non_zero_slider)
      zero_volcano_plot(data,
                   "Rank(Median)",
                  "Number of Zeros",
                  "Median Read Counts and Number of Zeros Over All Genes")}})
  
  output$count_heatmap <- renderPlot(
    {
      req(input$count_file)
      table <- load_count_data()
      if(!is.null(table)){
        data <- filter_count_data(table,
                                  input$variance_slider,
                                  input$non_zero_slider)
      
        
        t <- "Clustered Heatmap of Counts Remaining After Filtering"
        plot_heatmap(data, t)
        }}, height = 700, width = 900)
  
  output$count_pca <- renderPlot(
    {
      req(input$count_file)
      table <- load_count_data()
      if(!is.null(table)){
        data <- filter_count_data(table,
                                  input$variance_slider,
                                  input$non_zero_slider)
        
      if(input$pc_button == "PC vs PC" && input$PC1 != "..." && input$PC2 != "..."){
        plot_pc_v_pc(data, input$PC1, input$PC2)
      }
      else
      {
        plot_pca_beeswarm(data)
      }
    }})
  
  
  # observe({
  #   x <- input$pc_button
  #   
  #   # Can also set the label and select items
  #   if(x == "Item B"){
  #   updateRadioButtons(session, "pc_button2",
  #                      label = paste("radioButtons label", x),
  #                      choices = c("new choice", "second new choice")
  #   )}
  # })
  

  
}

# Run the application
shinyApp(ui = ui, server = server)
