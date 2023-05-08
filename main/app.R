#' This is the main RShiny application for this web app
#' Author: Justin Reimertz
#' 

# Set the maximum file upload size
options(shiny.maxRequestSize = 50 * 1024^2)
# Import functions and all dependent libraries
source("all_functions.R")

##### Define UI #####
ui <- fluidPage(
  # Set bootswatch theme
  theme = bslib::bs_theme(bootswatch = "pulse"),
  # Set up title and short purpose statement for the application
  titlePanel(
    h1("RNA-Seq Analysis", align = "center")),
  h3("Purpose", align = "center"),
  hr(),
  p("This is an app that can be used to perform various analyses on a given
     RNA-seq dataset. All that is required to use this app effectively is a 
    comma or tab delimited file with normalized normalized RNA-seq counts data
    and a comma or tab delimited file with sample information. Each of the pages
    in this app produces a different output and serves a different function.",
    align = "center"),
  br(),
  # This application involves a series of nested tabs to distinguish between 
  # pages with user inputs and the various functions that can be accomplished
  tabsetPanel(
    ##### First Page ####
    # The first page is exploration of the meta data associated with the RNA-seq
    # data provided to the app
    tabPanel("Sample Information Exploration",
             sidebarLayout(
               sidebarPanel(h4("Meta Data Exploration"),
                            p("On this page you can explore the sample information 
                              associated with the RNA-seq data"),
                            # Input accepts a comma or tab delimited file with
                            # sample information. The placeholder is set to the
                            # meta data file that this application was initially
                            # tested on, though any "tidy" meta data file can
                            # be inputted
                            fileInput(inputId = "metaFile", 
                                      label = "Load sample information file",
                                      accept = c("text/comma-separated-values",
                                                 "text/tab-separated-values",
                                                 "text/csv", "text/tsv",
                                                 ".csv",".tsv"),
                                      buttonLabel = "Browse", 
                                      placeholder = "GSE64810_metadata.txt"),
                            h4("Set Plot Parameters"),
                            p("Diagnostic plots can be generated for the 
                              uploaded sample information by choosing a column
                              to plot, a column to group the data by, and a
                              type of plot to make."),
                            # Prompts the user to select columns using the
                            # column names in the inputted meta data file to
                            # create plots with to visualize and explore the 
                            # data
                            selectInput(inputId = "column", 
                                        label = "Select column to plot",
                                        choices = "Pending Upload"),
                            selectInput(inputId = "group",
                                        label = "Select column to group by",
                                        choices = "Pending Upload"),
                            # Prompts the user to select a type of plot to 
                            # visualize the selected column and group with
                            radioButtons(inputId = "plotType", 
                                         label = "Choose the plot type",
                                         choices = c("Histogram","Density","Violin"), 
                                         selected = "Violin"),
                            # Plots will only be generated when this button is
                            # pressed
                            actionButton(inputId = "plotSpecs", 
                                         label = "Plot", 
                                         style = "color: white; background-color: #5b39a1",
                                         icon = icon("chart-column", 
                                                     class = "fa-solid fa-chart-column",
                                                     lib = "font-awesome"))
                            ),
               mainPanel(
                 # Each table displays a different output based on the inputted
                 # data and parameters set on the side panel
                 tabsetPanel(
                   # Outputs a summary table of the columns and data types 
                   # stored in the uploaded meta data file
                   tabPanel(title = "Data Summary", 
                            DT::dataTableOutput(outputId = "summaryTable")),
                   # Outputs the full meta data table
                   tabPanel(title = "Full Data Table", 
                            DT::dataTableOutput(outputId = "dataTable")),
                   # Outputs the selected plot type using the selected column
                   # and grouping variables
                   tabPanel(title = "Data Plots",
                            plotOutput(outputId = "metaPlot"))
                   )
                 )
               )
             ),
    ##### Second Page #####
    # The second page is an exploration of the normalized RNA-seq counts matrix
    # data file
    tabPanel("Counts Matrix Exploration",
             sidebarLayout(
              sidebarPanel(h4("Counts Data Exploration"),
                            p("On this page you can explore the counts matrix 
                              for the RNA-seq data. This will allow for 
                              filtering and visualization of the data."),
                           # Input accepts a comma or tab delimited file with
                           # sample information. The placeholder is set to the
                           # normalized counts data file that this application 
                           # was initially tested on, though any "tidy" 
                           # normalized counts data file can be inputted
                           fileInput(inputId = "countsFile", 
                                     label = "Load normalized counts file",
                                     accept = c("text/comma-separated-values",
                                                "text/tab-separated-values",
                                                "text/csv", "text/tsv",
                                                ".csv",".tsv"),
                                     buttonLabel = "Browse", 
                                     placeholder = "GSE64810_norm_counts_adjust.txt"),
                           tags$style(".btn-file {background-color: #5b39a1;}"),
                           # Prompts the user to select a variance threshold
                           # to filter genes by
                           h4("Data Filters"),
                           p("Use the sliders to set thresholds for filtering
                             the normalized count data"),
                           sliderInput(inputId = "varThreshold",
                                       label = "Percentile of variance between genes",
                                       min = 0, max = 1, value = 0.6),
                           # Prompts user to select a minimum number of samples
                           # with non-zero counts to filter genes by
                           sliderInput(inputId = "sampleThreshold",
                                       label = "Number of samples with non-zero counts",
                                       min = 0, max = 100, value = 60),
                           # Filters will only be applied to the counts data 
                           # when this button is pressed
                           actionButton(inputId = "filtersApplied", 
                                        label = "Apply Filters", 
                                        style = "color: white; background-color: #5b39a1",
                                        icon = icon("cubes", 
                                                    class = "fa-thin fa-cubes",
                                                    lib = "font-awesome")),
                           # Prompts the user to select the number of principal
                           # components to include in the PCA beeswarm plots
                           h4("Principal Component Analysis"),
                           p("Perform PCA on the normalized count data"),
                           # Only perform PCA after button is pressed
                           actionButton(inputId = "doPCA", 
                                        label = "Perform PCA", 
                                        style = "color: white; background-color: #5b39a1",
                                        icon = icon("down-left-and-up-right-to-center", 
                                                    class = "fa-sharp fa-light 
                                                    fa-down-left-and-up-right-to-center",
                                                    lib = "font-awesome")),
                           # After performing PCA output a slider to select the
                           # number of PCs to include in the beeswarm plot and
                           # a button to make the PCA plot
                           # Output principal component number slider
                           sliderInput(inputId = "numPC", 
                                         label = "Number of principal components 
                                                  to include in beeswarm plot", 
                                         min=0, max=3, value = 2),
                           actionButton(inputId = "makePCAPlot", 
                                          label = "Make PCA Plot", 
                                          style = "color: white; background-color: #5b39a1",
                                          icon = icon("chart-column", 
                                                      class = "fa-solid fa-chart-column",
                                                      lib = "font-awesome"))
                           ),
              mainPanel(
                # Each tab displays a different output based on the inputted
                # data and parameters set on the side panel
                tabsetPanel(
                  # Outputs a table summarizing the effects of the filters
                  # including the number and percentage of genes passing the
                  # current filters and the number and percentage of genes that
                  # were removed based on the applied filters
                  tabPanel(title = "Data Summary",
                           tableOutput(outputId = "filteredDataTable")
                           ),
                  # Outputs diagnostic scatter plots of the full data but
                  # distinguishing between genes that were filtered out and
                  # those that passed the current filters
                  tabPanel(title = "Diagnostic Scatter Plots",
                           fluidRow(column(2, align = "right", 
                                           plotOutput(outputId = "varPlot", 
                                                      width  = "600px",
                                                      height = "500px"),
                                           plotOutput(outputId = "samplePlot",
                                                      width  = "600px",
                                                      height = "500px")
                                           )
                                    )
                           ),
                  # Outputs a clustered heatmap of the genes passing the 
                  # current filters
                  tabPanel(title = "Clustered Heatmap",
                           plotOutput("filterHeatmap")),
                  # Outputs a beeswarm plot of the top N principal components
                  # selected by the user
                  tabPanel(title = "Principal Component Analysis",
                           plotOutput("PCAPlot"))
                  )
                )
              )
             ),
    ##### Third Page #####
    # The third page explores differential expression analysis results either
    # from a file uploaded by the user or performs differential expression 
    # analysis using DESeq2 on the uploaded normalized counts data
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(h4("Differential Expression"),
               p("On this page you can explore the results of differential 
               expression analysis. Either upload results of an already completed 
               analysis or differential expression can be performed using the 
               DESeq2 package."),
                            # If the user checks this box a fileupload input
                            # will appear so that the user can upload their 
                            # differential results file. This value is set to
                            # FALSE so that if the user does not check the box
                            # differential expression analysis will be performed
                            # on the upload counts data
                            checkboxInput(inputId = "diffExprs",
                                          label = "Upload differential 
                                          expression data?",
                                          value = FALSE),
                            uiOutput(outputId = "fileUpload"),
                            uiOutput(outputId = "pheno1Input"),
                            uiOutput(outputId = "pheno2Input"),
                            uiOutput(outputId = "diffExprsButton"),
                            h4("Volcano Plot"),
                            p("A volcano plot can be generated with",
                              span(glue('"log2 fold-change"'), style = "strong"),
                              "on the x-axis and",
                              span("p-adjusted", style = "strong"),
                              "on the y-axis"),
                            radioButtons(inputId = "x_axis", label = "Choose the 
                                         column for the x-axis",
                                         choices = c("baseMean","log2FoldChange",
                                                     "lfcSE","stat","pvalue",
                                                     "padj"), 
                                         selected = "log2FoldChange"),
                            radioButtons(inputId = "y_axis", label = "Choose the 
                                         column for the y-axis",
                                         choices = c("baseMean","log2FoldChange",
                                                     "lfcSE","stat","pvalue",
                                                     "padj"), 
                                         selected = "padj"),
                            colourInput(inputId = "base", 
                                        label = "Base point color", 
                                        value = "#771ec1"),
                            colourInput(inputId = "highlight", 
                                        label = "Highlight point color", 
                                        value = "#00e3db"),
                            sliderInput(inputId = "pValueMag",
                                        label = "Select the magnitude of the 
                                        adjusted p-value coloring",
                                        min = -50, max = 0, value = -10),
                            actionButton(inputId = "makeVolcano", label = "Plot", 
                                         style = "color: white; background-color: #5b39a1",
                                         icon = icon("chart-column", 
                                                     class = "fa-solid fa-chart-column",
                                                     lib = "font-awesome"))
                            ),
               mainPanel(
                 # Each tab displays a different output based on the inputted
                 # data and parameters set on the side panel
                 tabsetPanel(
                   # Outputs a table with the full differential expression
                   # analysis results
                   tabPanel(title = "Differential Expression Results",
                            dataTableOutput(outputId = "diffExprsResults")),
                   # Outputs a volcano plot with points colored by a significance
                   # threshold
                   tabPanel(title = "Volcano Plot",
                            plotOutput(outputId = "volcanoPlot"))
                   )
                 )
               )
             ),
    ##### Fourth Page #####
    # The fourth page performs pairwise correlation network analysis on an 
    # inputted list of genes found in the normalized counts data file and 
    # generates a clustered heatmap of these genes, as well as a network
    # correlation graph
    tabPanel("Correlation Network Analysis",
             sidebarLayout(
               sidebarPanel(h4("Correlation Network Analysis"),
                              p("On this page you can perform pairwise gene
                              expression correlation for a given set of genes.
                              In addition to inputting genes to perform this
                              analysis on, a minimum correlation threshold can
                              be set for an edge between two genes to be 
                              included."),
                            # Input accepts a list of HGNC gene names and/or 
                            # ensembl ids that are included in the normalized
                            # counts data file
                            textAreaInput(inputId = "geneList",
                                          label = "Input genes to inlcude in 
                                          pariwise correlation analysis",
                                          cols = 1,
                                          placeholder = "APOE"),
                            # Slider that sets the threshold for minimum 
                            # correlation for drawing an edge between two genes
                            sliderInput(inputId = "corrThreshold",
                                        label = "Minimum correlation threshold for
                                        drawing an edge",
                                        min = 0, max = 1, value = 0.3),
                            # Only perform pairwise correlation analysis 
                            # when this button is pressed
                            actionButton(inputId = "corrGenes", 
                                         label = "Perform Pairwise Correlation", 
                                         style = "color: white; background-color: #5b39a1",
                                         icon = icon("circle-nodes", 
                                                     class = "fa-light fa-circle-nodes",
                                                     lib = "font-awesome"))
                            ),
               mainPanel(
                 # Each tab displays a different output based on the inputted
                 # data and parameters set on the side panel
                 tabsetPanel(
                   # Outputs a clustered heatmap of the specified genes 
                   tabPanel(title = "Clustered Heatmap",
                            plotOutput(outputId = "corrHeatmap")),
                   # Outputs a correlation network graph
                   tabPanel(title = "Correlation Network",
                            plotOutput(outputId = "corrNetwork")),
                   # Outputs a table with the metrics determined by the 
                   # correlation network
                   tabPanel(title = "Correlation Metrics Table",
                            dataTableOutput(outputId = "corrResults"))
                   )
                 )
               )
             )
  )
)

##### Define server function #####
server <- function(input, output, session) {
  ##### Meta Data Output #####
  # Output meta data summary table
  output$summaryTable <- DT::renderDataTable({
    DT::datatable(make_summary_table(load_meta_data(input)()))
    })
  # Output the full meta data table
  output$dataTable <- DT::renderDataTable({
    DT::datatable(load_meta_data(input)())
  })
  # Update selectInputs based on column names of the inputted meta data file
  observe({
    updateSelectInput(session, "column", choices = names(load_meta_data(input)()))
  })
  observe({
    updateSelectInput(session, "group", choices = names(load_meta_data(input)()))
  })
  # Output the selected meta data plot
  output$metaPlot <- renderPlot({
    input$plotSpecs
    
    isolate({
      make_meta_plot(
        load_meta_data(input)(), input$column, input$group, input$plotType
      )})
  }, width = 800, height = 600)
  
  ##### Counts Data Output #####
  # Create variable storing the results from the currently set filters whenever
  # the 'Apply Filters' button is pressed
  filteredCounts <- eventReactive(input$filtersApplied,{
    get_filtered_counts(load_count_data(input)(), input$varThreshold,
                  input$sampleThreshold)
  })
  # Create variable for list of filtered genes
  filteredGenes <- eventReactive(input$filtersApplied,{
    filteredCounts()$gene
  })
  # Create variable for summary data after filters are applied
  filtersSummary <- eventReactive(input$filtersApplied,{
    apply_filters(load_count_data(input)(), input$varThreshold, 
                  input$sampleThreshold)
  })
  # Create a variable for the PCA results scores
  
  # Perform PCA and store results in a variable
  pcaResults <- eventReactive(input$doPCA,{
    # If filters have been applied to the count data use the filtered data
    if (!is.null(filteredCounts)) {
      do_PCA(filteredCounts()) 
    }
    # If no filters have been applied to the data use the full dataset
    else {
      do_PCA(load_count_data(input)())
    }
  })
  # Output counts data 
  # If filters have been applied output the text, update after 'Apply Filters'
  # button is pressed re-actively
  output$filteredDataTable <- renderTable({
    filtersSummary()
  })
  # Diagnostic scatter plots output
  # Output the median count vs variance plot when 'Apply Filters' button is 
  # pressed
  output$varPlot <- renderPlot({
    input$filtersApplied
    
    isolate({
      make_diag_plot(
        load_count_data(input)(), filteredCounts(), 
        filteredGenes(), "varianceScatter"
        )
      })
  }) 
  # Output the median count vs number of zeros plot when 'Apply Filters' button
  # is pressed
  output$samplePlot <- renderPlot({
    input$filtersApplied
    
    isolate({
      make_diag_plot(
        load_count_data(input)(), filteredCounts(),
        filteredGenes(), "sampleScatter"
      )
    })
  })
  # Output clustered heatmap
  # Output the clustered heatmap of filtered genes when 'Apply Filters' button
  # is pressed
  output$filterHeatmap <- renderPlot({
    input$filtersApplied
    
    isolate({
      make_heatmap(filteredCounts())
    })
  })
  observe({
    PCAlist <- pcaResults()
    updateSliderInput(session, "numPC", max = length(colnames(PCAlist$x)))
  })
  # Output PCA plot
  output$PCAPlot <- renderPlot({
    input$makePCAPlot
    
    isolate({
        PCA_plot(pcaResults(), input$numPC)
    })
  }, width = 1000, height = 800)
  ##### Differential Expression Output #####
  # Create global variable for differential expression results table if results
  # were uploaded
  deResultsTable <- eventReactive(input$showDEResults,{
    if (!input$diffExprs) {
      return(NULL)
    }
    else {
      load_de_results(input)() 
    }
  })
  # Create variable for DESeq2 results table
  DESeqResultsTable <- eventReactive(input$doDiffAnalysis,{
    if (input$diffExprs) {
      return(NULL)
    }
    else {
      run_deseq(filteredCounts(), load_meta_data(input)(), 
                c(input$pheno1, input$pheno2), 10)
    }
  })
  # If user checks box to upload previously computed differential expression
  # results render the file input 
  output$fileUpload <- renderUI({
    if (!input$diffExprs) {
      return(NULL)
    }
    fileInput(inputId = "deResultsFile", 
              label = "Load differential expression results file",
              accept = c("text/comma-separated-values",
                         "text/tab-separated-values",
                         "text/csv", "text/tsv",
                         ".csv",".tsv"),
              buttonLabel = "Browse", 
              placeholder = "Results File")
  })
  # Prompt the user to choose which phenotypes to perform the differential
  # expression comparison between
  output$pheno1Input <- renderUI({
    if (input$diffExprs) {
      return(NULL)
    }
    textInput(inputId = "pheno1", label = "Input a phenotype", value = "")
  })
  output$pheno2Input <- renderUI({
    if (input$diffExprs) {
      return(NULL)
    }
    textInput(inputId = "pheno2", label = "Input a phenotype", value = "")
  })
  # Change the button to reflect if the user is uploading differential
  # expression results or if differential expression should be calculated from
  # the uploaded counts data
  output$diffExprsButton <- renderUI({
    if (!input$diffExprs) {
      # Only perform differential expression analysis 
      # when this button is pressed
      actionButton(inputId = "doDiffAnalysis", 
                   label = "Do Differential Expression", 
                   style = "color: white; background-color: #5b39a1",
                   icon = icon("code-compare", 
                               class = "fa-light fa-code-compare",
                               lib = "font-awesome"))
    }
    else {
      # Only render the uploaded differential expression results when this 
      # button is pressed
      actionButton(inputId = "showDEResults", 
                   label = "Display Results", 
                   style = "color: white; background-color: #5b39a1",
                   icon = icon("code-compare", 
                               class = "fa-light fa-code-compare",
                               lib = "font-awesome"))
    }
  })
  # Render the data table of differential expression results
  output$diffExprsResults <- DT::renderDataTable({
    if (input$diffExprs) {
      DT::datatable(deResultsTable())
    }
    else {
      DT::datatable(DESeqResultsTable())
    }
  })
  # Display the volcano plot
  output$volcanoPlot <- renderPlot({
    input$makeVolcano
    
    isolate ({
      if (input$diffExprs) {
        volcano_plot(
          deResultsTable(), input$x_axis, input$y_axis, 
          input$pValueMag, input$base, input$highlight
        )
      }
      else {
        volcano_plot(
          DESeqResultsTable(), input$x_axis, input$y_axis, 
          input$pValueMag, input$base, input$highlight
        ) 
      }
    })
  }, width = 800, height = 600)
  ##### Correlation Network Output #####
  # Create variable with vector of inputted genes
  geneVector <- eventReactive(input$corrGenes,{
    unlist(strsplit(input$geneList, split = "[\r\n]"))
  })
  # Create a variable with ensembl gene ids and hgnc gene symbols stored in a
  # tibble for quick lookup
  geneLookup <- eventReactive(input$corrGenes,{
    geneTibble <- tibble(gene = geneVector())
    gene_match(load_count_data(input)(), geneTibble)
  })
  # Create a variable with the table of pairwise correlation metrics for the 
  # given list of genes
  corrTable <- eventReactive(input$corrGenes,{
    pairwise_correlation(geneLookup(), input$corrThreshold, "table")
  })
  # Perform pairwise correlation for the list of genes and output the result
  # as a clustered heatmap
  output$corrHeatmap <- renderPlot({
    input$corrGenes
    
    isolate({
      pairwise_correlation(geneLookup(), input$corrThreshold, "heatmap")
    })
  })
  # Perform pairwise correlation for the list of genes and output the resulting
  # correlation network graph
  output$corrNetwork <- renderPlot({
    input$corrGenes
    
    isolate({
      pairwise_correlation(geneLookup(), input$corrThreshold, "graph")
    })
  }, width = 800, height = 600)
  # Perform pairwise correlation for the list of genes and output a table with 
  # the resulting correlation metrics for each gene
  output$corrResults <- DT::renderDataTable({
    DT::datatable(corrTable())
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
