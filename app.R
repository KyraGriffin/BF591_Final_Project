## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7

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
    titlePanel("BF591 Assignment 7"),
    markdown(paste0("To use this application, download the CSV `deseq_res.csv`",
                    " from the data directory of this app's repository.")),
    
    # Sidebar with a slider input for p-value magnitude
    sidebarLayout(
        sidebarPanel(
            fileInput(
                "file",
                label = "Load differential expression results",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"),
                placeholder = "deseq_res.csv"
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
                value = "#22577A",
                closeOnClick = T
            ),
            colourInput(
                inputId = "highlight",
                label = "Highlight point color",
                value = "#FFCF56",
                closeOnClick = T
            ),
            sliderInput(
                "slider",
                "Select the magnitude of the p adjusted coloring:",
                min = -300,
                max = 0,
                value = -150,
                step = 1
            ),
            submitButton(
                text = "Plot",
                icon = icon("car-crash"),
                width = "100%"
            )
        ),
        # Show the volcano plot
        mainPanel(tabsetPanel(
            tabPanel("Plot", {
                plotOutput("volcano")
            }),
            tabPanel("Table",
                     tableOutput("table"))
        ))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
        df <- read.csv(input$file$datapath)
        colnames(df)[1] <- "gene"
        return(df)
    })
    
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
            p <- ggplot(dataf, aes(x = !!sym(x_name),
                           y = -log10(!!sym(y_name)))) +
                geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
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
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
        df_out <- dataf[which(dataf$padj < 1 * 10 ^ (as.numeric(slider))),]
        df_out$pvalue <- formatC(df_out$pvalue, digits = -2)
        df_out$padj <- formatC(df_out$padj, digits = -2)
        return(df_out)
    }
    
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    output$volcano <- renderPlot({
        req(input$file)
        df <- load_data()
        p <-volcano_plot(df,
                         input$x_axis,
                         input$y_axis,
                         input$slider,
                         input$base,
                         input$highlight)
        return(p)
    }, height = 700)
    
    # Same here, just return the table as you want to see it in the web page
    output$table <- renderTable({
        req(input$file)
        table <- load_data()
        colnames(table)[1] <- "gene"
        return(draw_table(dataf = table, slider = input$slider))
    }, striped = T)
}

# Run the application
shinyApp(ui = ui, server = server)
