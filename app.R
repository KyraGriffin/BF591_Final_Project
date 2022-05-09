## Author: Kyra Griffin
## ke31grif@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(DT)
library(ggbeeswarm)
library(plotly)
library(gplots)
library(qwraps2)
library(heatmaply)
library(shinythemes)
library(patchwork)
library(tidyverse)
library(bslib)
library(ggplot2)
library(plotly)
library(colourpicker)
# install.packages('shinyHeatmaply')
# install.packages("devtools")
# devtools::install_github("thomasp85/patchwork")
# install.packages("gplots")


# Radio Button Choices
dataset_choice <- c("Neurologically normal", "Huntington's Disease")

sample_y_choice <- c(
  "age_of_death", "AvgSpotLen", "Bases", "Bytes",
  "mrna.seq_reads", "pmi", "rin", "age_of_onset", "cag",
  "Duration", "h.v_cortical_score", "h.v_striatal_score",
  "vonsattel_grade")

deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")


# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("lumen"),
  tags$head(
    tags$style(HTML("
                  .btn {
                    height: 37px;
                    width: 100px;}"))
  ),
  titlePanel("BF591 Final Project - Kyra Griffin"),
  h4("Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"),
  tabsetPanel(
    tabPanel(
      "Samples",
      p("In this set of tabs we will be focusing on the metadata of the dataset. 
        On the side there is a file input that allows you to enter the sample metadata
        and choose which diagnosis you want to specifically look at."),
      br(),
      sidebarLayout(
        sidebarPanel(
          fileInput(
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
          submitButton(
            text = "Plot",
            icon = icon("car-crash"),
            width = "100%"
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Summary", 
                     p("The summary tab shows a summary table of the entire sample metadata with the column name, 
        type of the column, and either the mean(sd) or the distict values within a column.
        This table is also sorted but the diagnosis chosen."),br(),
                     tableOutput("samp_sum_table")),
            tabPanel("Table", 
                     p("The table tab shows a sortable and searchable data table of the entire 
                     sample metadata.This table is also sorted but the diagnosis chosen."),br(),
                     DT::dataTableOutput("sample_DT")),
            tabPanel("Plots", 
                     p("The plots tab shows violin plots of continuous variables within the sample metadata. 
                       These plots are not filtered by the diagnosis but when you chose age_of_onset,
                       cag,Duration, h.v_cortical_score, h.v_striatal_score, or vonsattel_grade only one
                       representation shows because these categories or values are only present for the 
                       patients with Huntington's Disease"),br(),
                     radioButtons(
              inputId = "samp_y_axis",
              inline = TRUE,
              label = "Choose the y-axis",
              choices = sample_y_choice,
              selected = "age_of_death"
            ), plotOutput("sample_plot"))
          ),
        )
      )
    ),
    tabPanel(
      "Counts",p("The Counts focuses on the normalized conts matrix. 
                 The side panel has a file upload area as well as two sliders 
                 for the user to chose the percentaile of variance and the 
                 number of non-zero samples a gene can contain."),br(),
      sidebarLayout(
        sidebarPanel(
          fileInput(
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
          submitButton(
            text = "Plot",
            icon = icon("car-crash"),
            width = "100%"
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Summary", br(),p("The summary tab shows a summary table of
                                       the number of samples total, the total number of genes, 
                                       the total number and the percentage of genes that pass 
                                       and that do not pass the filtering set by the sliders."),br(),
                     tableOutput("count_sum_table")),
            tabPanel(
              "Diagnostic Scatter Plot", br(),
              p("The Diagnostic Scatter Plot tab shows two plots,median count vs variance and 
              median count vs number of zeros, where genes passing filters are marked in a darker 
              color, and genes filtered out are lighter."),br(),
              splitLayout(
                cellWidths = c("50%", "50%"),
                plotOutput("median_vs_var_plot"),
                plotOutput("median_vs_zeros")
              )
            ),
            tabPanel(
              "Clustered Heatmap", br(),
              p("The Clustered Heatmap tab shows a clustered heatmap of counts remaining after 
                filtering, with a color bar located at the top left and the sample names on 
                the x-axis and the genes on the y-axis"),br(),
              plotOutput("count_heatmap")
            ),
            tabPanel(
              "PCA",
              p("The PCA tab shows either two principle components plotted against one 
                another or two Beeswarm plots of the top two PCs grouped by the diagnosis they represent."),br(),
              h4("Select PC vs PC or Beeswarm to plot either two principle components against
              one another or the top 2 principle components in Beeswarm Plots:"),
              radioButtons(
                "pc_button", "",
                c("PC vs PC", "Beeswarm Plot")
              ),
              h4("If you choose PC vs PC please enter the numbers of the two principle components you want to plot below:"),
              textInput("PC1", label = h6("First PC"), value = "...", width = 50),
              textInput("PC2", label = h6("Second PC"), value = "...", width = 70),
              plotOutput("count_pca")
            )
          )
        ),
      )
    ),
    tabPanel(
      "DE",
      p("The Differential Expression or DE tab focuses on the deseq results. 
        The side bar takes the deseq results as input and allows the user to 
        choose the x and y-axis, as well as the colors and the magnitude of the 
        adjusted p-value for the volcano plot that will be produced in the Plot tab."),br(),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "DE_file",
            label = "Load differential expression results",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            ),
            placeholder = "GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv"
          ),
          HTML(
            paste(
              "A volcano plot can be generated with",
              "<b>\"log<sub>2</sub> fold-change\"</b> on the x-axis and",
              "<b>\"p-adjusted\"</b> on the y-axis.<br>"
            )
          ),
          br(),
          radioButtons(
            inputId = "x_axis",
            label = "Choose the column for the x-axis",
            choices = deseq_choices,
            selected = "log2FoldChange"
          ),
          radioButtons(
            inputId = "y_axis",
            label = "Choose the column for the y-axis",
            choices = deseq_choices,
            selected = "padj"
          ),
          colourInput(
            inputId = "base",
            label = "Base point color",
            value = "#F8766D",
            closeOnClick = T
          ),
          colourInput(
            inputId = "highlight",
            label = "Highlight point color",
            value = "#619CFF",
            closeOnClick = T
          ),
          sliderInput(
            "slider",
            "Select the magnitude of the p adjusted coloring:",
            min = -35,
            max = 0,
            value = -10,
            step = 1
          ),
          submitButton(
            text = "Plot",
            icon = icon("car-crash"),
            width = "100%"
          )
        ),
        # Show the volcano plot
        mainPanel(
          tabsetPanel(
            tabPanel("Data Table",
                     p("The Data Table tab shows sortable and searchable table displaying 
                       the differential expression results."),br(),
                     DT::dataTableOutput("DE_DT")),
            tabPanel("Plot", 
                     p("The Plot tab shows a plot of the DE results with the points 
                       where the padj is less than 1 * 10^ magnitude set in the side bar colored in the highlight color selected."),br(),
                     plotOutput("volcano")),
            tabPanel("Plot Related Table",
                     p("The Plot Related Table tab shows a table of the genes 
                       that have a padj is less than 1 * 10^ magnitude set, corresponding to the results of the Plot tab."),br(),
                     tableOutput("DE_plot_table"))
          )
        ),
      )
    ),
    tabPanel(
      "GSEA",
      p("The GSEA tab focuses on analyzing or looking closer at the fgsea results 
        generated in run_fgsea.R. This set of tabs has a side bar that allows for 
        data input of the fgsea.csv results file generated."),br(),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "GSEA_file",
            label = "Load FGSEA file:",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            ),
            placeholder = "fgsea.csv"
          ),
          submitButton(
            text = "Submit",
            icon = icon("car-crash"),
            width = "100%"
          )
        ),
        # Show the volcano plot
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Top Results",
              p("The Top Results tab takes a slider input of the number of top pathways
                that the user wants plotted. The bar plot of the results is 
                displayed to the right of the slider input."),br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    "GSEA_slider",
                    "Select the number of top pathways to plot by adjusted p-value:",
                    min = 0,
                    max = 30,
                    value = 15,
                    step = 1
                  ),
                  submitButton(
                    text = "Submit",
                    icon = icon("car-crash"),
                    width = "100%"
                  )
                ),
                mainPanel(
                  plotOutput("top_paths_plot")
                ),
              )
            ),
            tabPanel(
              "Table",
              p("The Table tab takes a slider input of adjusted pvalue threshold.
                As well as whether the user wants to see only the posotive, 
                negative, or all of the NES pathways. The table is filtered based
                on the slider and pathway selection and is sortable and searchable"),br(),
                p("This filtered table can also be dowloaded into a CSV file using the download button under
                  the submit button."),br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    inputId = "GSEA_p_thresh_slider",
                    label = "Select the adjusted p-value threshold:",
                    min = 0,
                    max = 1,
                    value = 0.1,
                    step = 0.1
                  ),
                  radioButtons(
                    inputId = "pathways",
                    label = "Select whether you want Positive, Negative, or All NES pathways:",
                    choices = c("Positive", "Negative", "All"),
                    selected = "All"
                  ),
                  submitButton(
                    text = "Submit",
                    icon = icon("car-crash"),
                    width = "100%"
                  ),
                  br(),
                  downloadButton("downloadData", "Download ", width = "100%")
                ),
                mainPanel(
                  DT::dataTableOutput("gsea_table")
                ),
              )
            ),
            tabPanel(
              "FGSEA Plot",
              p("The FGSEA Plot tab takes a slider input to filter table by 
              adjusted p-value. A scatter plot of NES on x-axis and -log10 
              adjusted p-value on y-axis, with gene sets below threshold in grey color 
                is generated and reactive to the slider input."),br(),
              sidebarLayout(
                sidebarPanel(
                  sliderInput(
                    inputId = "GSEA_thresh_plot_slider",
                    label = "Select the adjusted p-value threshold:",
                    min = -30,
                    max = 5,
                    value = -10,
                    step = 1
                  ),
                  submitButton(
                    text = "Submit",
                    icon = icon("car-crash"),
                    width = "100%"
                  )
                ),
                mainPanel(
                  plotOutput("gsea_plot")
                ),
              )
            ),
            tabPanel(
              "FGSEA Plot Table",
              p("The FGSEA Plot Table tab shows a table of the points or 
                pathways/gene sets below the threshold (the points show in grey) 
                set in the FGSEA Plot tab."),br(),
                mainPanel(
                  tableOutput("gsea_plot_table")
                ),
              
            ),
          )
        ),
      ),
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 35 * 1024^2)

  ############### Load Data ####################

  #' Load Sample Data
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
    df <- read.csv(file = inFile$datapath, header = TRUE, stringsAsFactors = TRUE)

    return(df)
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

  #' Load DE Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is
  #' still a "function", but it will take no arguments. The `reactive({})` bit
  #' says "if any of my inputs (as in, input$...) are changed, run me again".
  #' This is useful when a user clicks a new button or loads a new file. In
  #' our case, look for the uploaded file's datapath argument and load it with
  #' read.csv. Return this data frame in the normal return() style.
  load_DE_data <- reactive({
    df <- read_tsv(input$DE_file$datapath, col_names = TRUE, show_col_types = FALSE)
    colnames(df)[1] <- "gene"
    return(df)
  })
  
  #' Load FGSEA Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is
  #' still a "function", but it will take no arguments. The `reactive({})` bit
  #' says "if any of my inputs (as in, input$...) are changed, run me again".
  #' This is useful when a user clicks a new button or loads a new file. In
  #' our case, look for the uploaded file's datapath argument and load it with
  #' read.csv. Return this data frame in the normal return() style.
  load_GSEA_data <- reactive({
    inFile <- input$GSEA_file

    if (is.null(inFile)) {
      return(NULL)
    }
    df <- read.csv(file = inFile$datapath, header = TRUE, stringsAsFactors = TRUE)

    return(df)
  })

  ############### Functions ####################


  #' Draw summary table
  #'
  #' @param dataf Data frame loaded by load_sample_data()
  #' @param diagnosis Diagnosis of sample row from the radio button input.
  #'
  #' @return Data frame
  #'
  #' @examples draw_sum_table(sample_metadata, "Neurologically normal")
  draw_sum_table <- function(dataf, diagnosis) {
    df_out <- dataf %>%
      dplyr::filter(Diagnosis == diagnosis)

    if (diagnosis == "Neurologically normal") {
      df_out <- dataf %>%
        dplyr::select(-c(
          age_of_onset, cag, Duration, h.v_cortical_score,
          h.v_striatal_score, vonsattel_grade
        ))
    }
    # Manually looping through the columns of the filtered data 
    # to generate the vectors for the summary table.
    
    col_type <- c()
    dv_ms <- c()

    for (i in 1:ncol(df_out)) {
      col <- df_out[, i]
      type <- class(col)
      col_type <- append(col_type, type)

      if (type == "factor") {
        u_vals <- as.character(unique(df_out[, i]))
        dv_ms <- append(dv_ms, paste(u_vals, collapse = ", "))
      } else {
        m <- mean(df_out[, i])
        s <- sd(df_out[, i])
        mean_sd <- as.character(paste0(m, " (+-", s, ")"))
        dv_ms <- append(dv_ms, mean_sd)
      }
    }
    
    # Forming the data frame
    column_names <- colnames(df_out)
    type_col <- col_type
    mean_sd_col <- dv_ms
    df <- data.frame(column_names, type_col, mean_sd_col)
    colnames(df) <- c("Column Name", "Type", "Mean (sd) or Distinct Values")

    return(df)
  }

  #' Draw sample plot
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param y_val The category to represent the y axis.
  #'
  #' @return Data frame
  #'
  #' @examples draw_sample_plot(sample_metadata, "age_of_death")
  draw_sample_plot <- function(dataf, y_val) {
    df_out <- dataf
    hd <- c("age_of_onset", "cag", "Duration", "h.v_cortical_score", "h.v_striatal_score", "vonsattel_grade")

    # Selecting numeric or continuous categories
    plot_data <- df_out %>% dplyr::select(where(is.numeric))
    plot_data <- plot_data %>%
      add_column(Diagnosis = df_out$Diagnosis)
    
    # Filtering data by diagnosis if categories related to Huntington's Disease are chosen
    if (y_val %in% hd) {
      plot_data <- plot_data %>% dplyr::filter(Diagnosis == "Huntington's Disease")
    }
  
    # Generation violin plot of column
    p <- ggplot(plot_data, aes(x = Diagnosis, y = !!sym(y_val), fill = Diagnosis)) +
      geom_violin()
    theme_bw() +
      theme(legend.position = "bottom")

    return(p)
  }

  #' Draw count summary table
  #'
  #' @param count_data Data frame loaded by load_count__data()
  #' @param var_val Variance percentile user chose
  #' @param zero_val Number of samples that are non-zero per gene
  #'
  #' @return Data frame
  #'
  #' @examples draw_count_sum_table(norm_counts, 0.20, 10)
  draw_count_sum_table <- function(count_data, var_val, zero_val) {
    # number of samples
    num_samples <- NCOL(count_data)
    # total number of genes
    total_genes <- length(count_data$gene)

    # Filter data based on variance percentile and nonzeros
    nonzero_genes <- rowSums(count_data[-1] == 0) >= zero_val
    nonzero_counts <- count_data[nonzero_genes, ]
    zero_counts <- filter(count_data, !gene %in% nonzero_counts$gene)
    nonzero_counts$variance <- apply(nonzero_counts[, -c(1)], 1, var)

    # Calculation percentile
    percent <- quantile(nonzero_counts$variance, prob = var_val / 100)

    var_counts <- dplyr::filter(nonzero_counts, variance >= percent) %>%
      dplyr::select(-variance)

    # number and % of genes passing current filter
    num_pass <- length(var_counts$gene)
    n1 <- round((num_pass / total_genes) * 100, 2)
    pass <- c(num_pass, n1)
    num_pass <- toString(pass, width = 30)

    # number and % of genes not passing current filter
    np <- length(zero_counts$gene) + (length(nonzero_counts$gene) - length(var_counts$gene))
    n2 <- round(((np / total_genes) * 100), 2)
    not_pass <- c(np, n2)

    num_not_pass <- toString(not_pass, width = 30)

    df <- data.frame(num_samples, total_genes, num_pass, num_not_pass)
    colnames(df) <- c(
      "Number of Samples", "Total Number of Genes",
      "Number , % of genes passing the filter",
      "Number , % of genes not passing the filter"
    )
    return(df)
  }


  #' Filter Count data
  #'
  #' @param count_data Data frame loaded by load_count__data()
  #' @param var_val Variance percentile user chose
  #' @param zero_val Number of samples that are non-zero per gene
  #'
  #' @return Data frame`
  #'
  #' @examples filter_count_data(norm_counts, 0.20, 10)
  filter_count_data <- function(count_data, var_val, zero_val) {

    # Filter data based on variance percentile and nonzeros
    nonzero_genes <- rowSums(count_data[-1] == 0) >= zero_val
    nonzero_counts <- count_data[nonzero_genes, ]

    zero_counts <- filter(count_data, !gene %in% nonzero_counts$gene)
    nonzero_counts$variance <- apply(nonzero_counts[, -c(1)], 1, var)

    percent <- quantile(nonzero_counts$variance, prob = var_val / 100)
    var_counts <- dplyr::filter(nonzero_counts, variance >= percent) %>%
      dplyr::select(-variance)
    not_var_counts <- filter(nonzero_counts, !gene %in% var_counts$gene) %>%
      dplyr::select(-variance)
  
    # Labeling genes that pass and do not pass filter for easier plotting
    not_filtered <- rbind(zero_counts, not_var_counts) %>% mutate(volcano = "Not Filtered")
    filtered <- var_counts %>% mutate(volcano = "Filtered")

    data <- as.data.frame(rbind(filtered, not_filtered))


    return(data)
  }

  #' median count vs variance diagnostic scatter plot
  #'
  #' @param data Data frame loaded by load_count__data() and filtered by filter_count_data()
  #' @param x_lab String label for x-axis
  #' @param y_lab String label for y-axis
  #' @param title String label for title
  #'
  #' @return ggplot object
  #'
  #' @examples var_volcano_plot(filtered_count_data, "Rank(Median)","Variance","Median Read Counts and Variability of Zeros Over All Genes")
  var_volcano_plot <- function(data, x_lab, y_lab, title) {
    dt <- data %>% dplyr::select(-volcano)
    medians <- apply(dt[, -c(1)], 1, median)
    variances <- apply(dt[, -c(1)], 1, var)

    plot_data <- tibble::tibble(median = medians, variance = variances, volcano = data$volcano)
    plot_data$rank <- rank(plot_data$median)

    p <- ggplot(plot_data, aes(x = rank, y = variance)) +
      geom_point(aes(color = volcano)) +
      scale_color_manual(values = c("#619CFF", "#F8766D")) +
      theme_bw() +
      ggplot2::xlab(x_lab) +
      ggplot2::ylab(y_lab) +
      ggtitle(title)
    theme(legend.position = "bottom")

    return(p)
  }
  
  #' median count vs number of zeros diagnostic scatter plot
  #'
  #' @param data Data frame loaded by load_count__data() and filtered by filter_count_data()
  #' @param x_lab String label for x-axis
  #' @param y_lab String label for y-axis
  #' @param title String label for title
  #'
  #' @return ggplot object
  #'
  #' @examples zero_volcano_plot(filtered_count_data, "Rank(Median)","Number of Zeros", "Median Read Counts and Number of Zeros Over All Genes"
  zero_volcano_plot <- function(data, x_lab, y_lab, title) {
    dt <- data %>% dplyr::select(-volcano)
    medians <- apply(dt[, -c(1)], 1, median)
    num_zero <- rowSums(dt[, -c(1)] == 0)

    plot_data <- tibble::tibble(median = medians, zeros = num_zero, volcano = data$volcano)
    plot_data$rank <- rank(plot_data$median)

    p <- ggplot(plot_data, aes(x = rank, y = zeros)) +
      geom_point(aes(color = volcano)) +
      scale_color_manual(values = c("#619CFF", "#F8766D")) +
      theme_bw() +
      ggplot2::xlab(x_lab) +
      ggplot2::ylab(y_lab) +
      ggtitle(title)
    theme(legend.position = "bottom")

    return(p)
  }

  
  #' Clustered Heatmap
  #'
  #' @param data_filtered Data frame loaded by load_count__data() and filtered by filter_count_data()
  #' @param title String label for title
  #'
  #' @return heatmap object
  #'
  #' @examples plot_heatmap(filtered_count_data, "Counts Heatmap")
  plot_heatmap <- function(data_filtered, title) {
    d_filter <- data_filtered %>%
      dplyr::filter(volcano != "Not Filtered") %>%
      dplyr::select(-volcano) %>%
      column_to_rownames(var = "gene")

    # print(select_if(d_filter, is.numeric))

    p <- heatmap.2(as.matrix(d_filter),
      margins = c(5, 5),
      main = title,
      cexRow = 0.7,
      cexCol = 0.7, scale = "row", trace = "none"
    )


    return(p)
  }

  #' PC vs PC plot
  #'
  #' @param filtered_data Data frame loaded by load_count__data() and filtered by filter_count_data()
  #' @param metadata Sample metadata
  #' @param PC_1 String of a number represnting a PC
  #' @param PC_2 String of a number represnting a PC
  #'
  #' @return ggplot object
  #'
  #' @examples plot_pc_v_pc(filtered_count_data, metadata, "1", "2")
  plot_pc_v_pc <- function(filtered_data, metadata, PC_1, PC_2) {
    filtered_data <- filtered_data %>% # dplyr::filter(volcano == "Filtered") %>%
      dplyr::select(-volcano) %>%
      column_to_rownames(var = "gene")

    #print(select_if(filtered_data, is.numeric))

    pca <- prcomp(t(filtered_data))
    plot_data <- metadata

    pc_1 <- as.numeric(PC_1)
    pc_2 <- as.numeric(PC_2)


    plot_data$PC1 <- pca$x[, pc_1]
    plot_data$PC2 <- pca$x[, pc_2]
    percent_var <- pca$sdev^2 / sum(pca$sdev^2)

    title <- (paste0(
      "PC", pc_1, ": ", round(percent_var[pc_1] * 100), "% variance", " Vs. ",
      "PC", pc_2, ": ", round(percent_var[pc_2] * 100), "% variance"
    ))

    pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, col = Diagnosis)) +
      geom_point() +
      xlab(paste0("PC", pc_1, ": ", round(percent_var[pc_1] * 100), "% variance ")) +
      ylab(paste0("PC", pc_2, ": ", round(percent_var[pc_2] * 100), "% variance")) +
      ggtitle(title)

    return(pca_plot)
  }

  #' Beeswarm Plot PCA
  #'
  #' @param filtered_data Data frame loaded by load_count__data() and filtered by filter_count_data()
  #' @param metadata Sample metadata
  #'
  #' @return ggplot object
  #'
  #' @examples plot_pca_beeswarm(filtered_count_data, metadata)
  plot_pca_beeswarm <- function(filtered_data, metadata) {
    filtered_data <- filtered_data %>% # dplyr::filter(volcano == "Filtered") %>%
      dplyr::select(-volcano) %>%
      column_to_rownames(var = "gene")

    #print(select_if(filtered_data, is.numeric))

    pca <- prcomp(t(filtered_data))
    plot_data <- metadata

    pc_1 <- as.numeric("1")
    pc_2 <- as.numeric("2")


    plot_data$PC1 <- pca$x[, pc_1]
    plot_data$PC2 <- pca$x[, pc_2]
    percent_var <- pca$sdev^2 / sum(pca$sdev^2)


    pca_plot1 <- ggplot(plot_data) +
      geom_beeswarm(aes(x = Diagnosis, y = PC1, color = Diagnosis), cex = 2, size = 2, show.legend = FALSE) +
      ylab(paste0("PC", pc_1, ": ", round(percent_var[pc_1] * 100), "% variance")) +
      ggtitle("Beeswarm Plot of PC 1")

    pca_plot2 <- ggplot(plot_data) +
      geom_beeswarm(aes(x = Diagnosis, y = PC2, color = Diagnosis), cex = 2, size = 2) +
      ylab(paste0("PC", pc_2, ": ", round(percent_var[pc_2] * 100), "% variance")) +
      ggtitle("Beeswarm Plot of PC 2")

    count_pca_plot <- pca_plot1 | pca_plot2

    return(count_pca_plot)
  }


  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  #' @details I bet you're tired of these plots by now. Me too, don't worry.
  #' This is _just_ a normal function. No reactivity, no bells, no whistles.
  #' Write a normal volcano plot using geom_point, and integrate all the above
  #' values into it as shown in the example app. The testing script will treat
  #' this as a normal function.
  #'
  #' !!sym() may be required to access column names in ggplot aes().
  #'
  #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      p <- ggplot(dataf, aes(
        x = !!sym(x_name),
        y = -log10(!!sym(y_name))
      )) +
        geom_point(aes(color = !!sym(y_name) < 1 * 10^(as.numeric(slider)))) +
        theme_bw() +
        scale_color_manual(values = c(color1, color2)) +
        theme(legend.position = "bottom") +
        labs(color = paste0(y_name, " < 1 × 10^", slider))
      return(p)
    }

  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will
  #' evaluate it normally. Not only does this function filter the data frame to
  #' rows that are above the slider magnitude, it should also change the format
  #' of the p-value columns to display more digits. This is so that it looks
  #' better when displayed on the web page. I would suggest the function
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #' gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #' gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_DE_table <- function(dataf, slider) {
    df_out <- dataf[which(dataf$padj < 1 * 10^(as.numeric(slider))), ]
    df_out$pvalue <- formatC(df_out$pvalue, digits = -2)
    df_out$padj <- formatC(df_out$padj, digits = -2)
    return(df_out)
  }


  #' Function to plot top ten positive NES and top ten negative NES pathways
  #' in a barchart
  #'
  #' @param fgsea_results (tibble): the fgsea results in tibble format returned by
  #'   the previous function
  #' @param num_paths (int): the number of pathways for each direction (top or
  #'   down) to include in the plot. Set this at 10.
  #'
  #' @return ggplot with a barchart showing the top twenty pathways ranked by positive
  #' and negative NES
  #' @export
  #'
  #' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
  top_pathways <- function(fgsea_results, num_paths) {
    num_paths <- as.numeric(num_paths)

    top_pos <- fgsea_results %>%
      slice_max(NES, n = num_paths) %>%
      pull(pathway)
    top_neg <- fgsea_results %>%
      slice_min(NES, n = num_paths) %>%
      pull(pathway)

    subset <- fgsea_results %>%
      filter(pathway %in% c(top_pos, top_neg)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(plot_name = str_replace_all(pathway, "_", " "))

    title <- (paste0("Top ", num_paths, " Up-regulated and Down-regulated Pathways - FGSEA Results "))

    plot <- subset %>%
      mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
      ggplot() +
      geom_bar(aes(x = plot_name, y = NES, fill = NES > 0), stat = "identity", show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = "#F8766D", "FALSE" = "#619CFF")) +
      # theme_minimal(base_size = 8) +
      ggtitle(title) +
      theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 16)
      ) +
      ylab("Normalized Enrichment Score (NES)") +
      xlab("") +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
      coord_flip()

    return(plot)
  }

  #' Function to filter NES pathways based on padj threshold and pathway selector.
  #'
  #' @param gsea_data (tibble): the fgsea results in tibble format returned by
  #'   the previous function
  #' @param p_thresh: padj threshold
  #' @param pathways: String represnting type of pathways.
  #'
  #' @return data frame
  #'
  #' @examples fgsea_filtered_table <- gsea_filter_table(gsea_data, 0.1, "Positive")
  gsea_filter_table <- function(gsea_data, p_thresh, pathways) {
    gsea_filtered <- gsea_data %>%
      dplyr::mutate(status = case_when(
        NES > 0 ~ "UP",
        NES < 0 ~ "DOWN",
        TRUE ~ "NS"
      )) %>%
      dplyr::filter(padj < p_thresh)
    if (pathways == "Positive") {
      gsea_filtered <- gsea_filtered %>%
        dplyr::filter(status == "UP") %>%
        dplyr::select(-status)
    } else if (pathways == "Negative") {
      gsea_filtered <- gsea_filtered %>%
        dplyr::filter(status == "DOWN") %>%
        dplyr::select(-status)
    } else if (pathways == "All") {
      gsea_filtered <- gsea_filtered %>%
        dplyr::select(-status)
    }

    return(gsea_filtered)
  }

  #' Function to filter NES pathways based on padj threshold and pathway selector.
  #'
  #' @param gsea_data (tibble): the fgsea results in tibble format returned by
  #'   the previous function
  #' @param p_thresh: padj threshold
  #'
  #' @return ggplot object
  #'
  #' @examples gsea_plot(gsea_data, -10)
  gsea_plot <- function(gsea_data, p_thresh) {
    p <- ggplot(gsea_data, aes(
      x = NES,
      y = -log10(padj)
    )) +
      geom_point(aes(color = padj < 1 * 10^(as.numeric(p_thresh)))) +
      theme_bw() +
      scale_color_manual(values = c("#F8766D", "grey")) +
      ggtitle("NES Vs. -log10 adjusted p-value") +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 20)
      ) +
      labs(color = paste0("padj < 1 × 10^", (as.numeric(p_thresh))))

    return(p)
  }

  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will
  #' evaluate it normally. Not only does this function filter the data frame to
  #' rows that are above the slider magnitude, it should also change the format
  #' of the p-value columns to display more digits. This is so that it looks
  #' better when displayed on the web page. I would suggest the function
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #' gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #' gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_gsea_plot_table <- function(dataf, slider) {
    df_out <- dataf[which(dataf$padj < 1 * 10^(as.numeric(slider))), ]
    
    df_out$pval <- formatC(df_out$pval, digits = -2)
    df_out$padj <- formatC(df_out$padj, digits = -2)
    return(df_out)
  }


  ############### Output ####################
  # The following code generates the output seen in the ui by calling the above functions

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
    if (input$diagnosis_samp == "Neurologically normal") {
      table <- table %>%
        dplyr::select(-c(age_of_onset, cag, Duration, h.v_cortical_score, h.v_striatal_score, vonsattel_grade))
    }
    DT::datatable(table, options = list(orderClasses = TRUE))
  })

  output$sample_plot <- renderPlot({
    req(input$sample_metadata)
    table <- load_sample_data()
    if (!is.null(table)) {
      draw_sample_plot(
        table,
        input$samp_y_axis
      )
    }
  })

  output$count_sum_table <- renderTable(
    {
      req(input$count_file)
      table <- load_count_data()
      if (!is.null(table)) {
        draw_count_sum_table(
          table,
          input$variance_slider,
          input$non_zero_slider
        )
      }
    },
    spacing = "m",
    bordered = T
  )

  output$median_vs_var_plot <- renderPlot({
    req(input$count_file)
    table <- load_count_data()
    if (!is.null(table)) {
      data <- filter_count_data(
        table,
        input$variance_slider,
        input$non_zero_slider
      )
      var_volcano_plot(
        data,
        "Rank(Median)",
        "Variance",
        "Median Read Counts and Variability Over All Genes"
      )
    }
  })


  output$median_vs_zeros <- renderPlot({
    req(input$count_file)
    table <- load_count_data()
    if (!is.null(table)) {
      data <- filter_count_data(
        table,
        input$variance_slider,
        input$non_zero_slider
      )
      zero_volcano_plot(
        data,
        "Rank(Median)",
        "Number of Zeros",
        "Median Read Counts and Number of Zeros Over All Genes"
      )
    }
  })

  output$count_heatmap <- renderPlot(
    {
      req(input$count_file)
      table <- load_count_data()
      if (!is.null(table)) {
        data <- filter_count_data(
          table,
          input$variance_slider,
          input$non_zero_slider
        )


        t <- "Clustered Heatmap of Counts Remaining After Filtering"
        plot_heatmap(data, t)
      }
    },
    height = 700,
    width = 900
  )

  output$count_pca <- renderPlot({
    req(input$count_file)
    table <- load_count_data()
    if (!is.null(table)) {
      data <- filter_count_data(
        table,
        input$variance_slider,
        input$non_zero_slider
      )

      if (input$pc_button == "PC vs PC" && input$PC1 != "..." && input$PC2 != "...") {
        metadata <- read.csv(file = "data/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
        plot_pc_v_pc(data, metadata, input$PC1, input$PC2)
      } else if (input$pc_button == "Beeswarm Plot") {
        metadata <- read.csv(file = "data/sample_metadata.csv", header = TRUE, stringsAsFactors = TRUE)
        plot_pca_beeswarm(data, metadata)
      }
    }
  })

  output$DE_DT <- DT::renderDataTable({
    req(input$DE_file)
    table <- load_DE_data()

    DT::datatable(table, options = list(orderClasses = TRUE))
  })

  output$volcano <- renderPlot(
    {
      req(input$DE_file)
      df <- load_DE_data()
      p <- volcano_plot(
        df,
        input$x_axis,
        input$y_axis,
        input$slider,
        input$base,
        input$highlight
      )
      return(p)
    },
    height = 700
  )

  output$DE_plot_table <- renderTable(
    {
      req(input$DE_file)
      table <- load_DE_data()
      colnames(table)[1] <- "gene"
      return(draw_DE_table(dataf = table, slider = input$slider))
    },
    striped = T
  )

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("dataset-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$GSEA_file)
      table <- load_GSEA_data()
      filtered_data <- gsea_filter_table(table, input$GSEA_p_thresh_slider, input$pathways)

      write.csv(filtered_data, file)
    }
  )

  output$top_paths_plot <- renderPlot(
    {
      req(input$GSEA_file)
      table <- load_GSEA_data()
      top_pathways(table, input$GSEA_slider)
    },
    height = 700,
    width = 900
  )

  output$gsea_table <- DT::renderDataTable(
    {
      req(input$GSEA_file)
      table <- load_GSEA_data()
      filtered_data <- gsea_filter_table(table, input$GSEA_p_thresh_slider, input$pathways)

      DT::datatable(filtered_data,
        extensions = "Buttons",
        options = list(
          orderClasses = TRUE
        )
      )
    },
    height = 700,
    width = 900
  )

  output$gsea_plot <- renderPlot(
    {
      req(input$GSEA_file)
      gsea_data <- load_GSEA_data()
      gsea_plot(gsea_data, input$GSEA_thresh_plot_slider)
    },
    height = 700,
    width = 900
  )

  output$gsea_plot_table <- renderTable(
    {
      req(input$GSEA_file)
      table <- load_GSEA_data()

      draw_gsea_plot_table(table, input$GSEA_thresh_plot_slider)
    }, striped = T,
    bordered = T
  )
  
  
}

# Run the application
shinyApp(ui = ui, server = server)