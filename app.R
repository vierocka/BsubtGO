library(shiny)
library(DT)
library(topGO)
library(clusterProfiler)

# Set working directory if necessary
setwd("~/Dropbox/Melih_AGMeier") # Modify the path to your files

# Load data outside the server function to improve performance
allGOwithConv <- read.table("PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv", sep="\t", header = TRUE)

# Prepare the gene universe and gene-to-GO mappings once
geneUniverse <- allGOwithConv$IDbsub

# Prepare gene-to-GO mappings for each ontology
prepareGene2GO <- function(data, column_index) {
  geneID2GOdata <- data[, c(2, column_index)]
  
  # Function to split GO terms and associate with gene IDs
  myExtensionFun <- function(x) {
    GOlist <- unlist(strsplit(as.character(x[2]), split = " "))
    IDarray <- rep(x[1], length(GOlist))
    ID_GO <- cbind(IDarray, GOlist)
    return(ID_GO)
  }
  
  # Apply the function to each row
  result_list <- apply(geneID2GOdata, 1, myExtensionFun)
  
  # Combine the results into a data frame
  ID_GO_df <- do.call(rbind, result_list)
  ID_GO_df <- as.data.frame(ID_GO_df)
  colnames(ID_GO_df) <- c("GeneID", "GOterm")
  
  # Create gene-to-GO mapping
  geneID2GO <- by(ID_GO_df$GOterm, ID_GO_df$GeneID, function(x) as.character(x))
  
  # Ensure that the names of geneID2GO match your gene universe
  geneID2GO <- geneID2GO[names(geneID2GO) %in% geneUniverse]
  
  return(geneID2GO)
}

# Prepare mappings for BP, MF, and CC
geneID2GObp <- prepareGene2GO(allGOwithConv, 3)
geneID2GOmf <- prepareGene2GO(allGOwithConv, 4)
geneID2GOcc <- prepareGene2GO(allGOwithConv, 5)

# Define UI for the application
ui <- fluidPage(
  # Application title
  titlePanel("Gene Ontology (GO) Term Analysis of Bacillus subtilis"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        inputId = "id_list",
        label = "Enter Gene IDs separated by commas:",
        placeholder = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300",
        rows = 6
      ),
      selectInput(
        inputId = "keytype",
        label = "Select strain of Bacillus subtilis:",
        choices = c("PG10", "BSUB168"),
        selected = "PG10"
      ),
      actionButton("analyze", "Analyze", class = "btn-primary"),
      br(),
      br(),
      downloadButton("downloadData", "Download Results")
    ),
    
    mainPanel(
      h4("Biological Processes:"),
      br(),
      DTOutput("result_tableBP"),
      br(),
      br(),
      h4("Molecular Function:"),
      br(),
      DTOutput("result_tableMF"),
      br(),
      br(),
      h4("Cellular Components:"),
      br(),
      DTOutput("result_tableCC")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to process input IDs and perform analysis
  analysis_results <- eventReactive(input$analyze, {
    req(input$id_list)  # Ensure input is available
    
    # Split the input IDs into a vector
    id_vector <- unlist(strsplit(input$id_list, split = ",|\\s+"))
    id_vector <- trimws(id_vector)  # Remove leading/trailing whitespace
    id_vector <- id_vector[nzchar(id_vector)]  # Remove empty strings
    
    validate(
      need(length(id_vector) > 0, "Please enter at least one gene ID.")
    )
    
    # Ensure IDs are in the gene universe
    genesOfInterest <- id_vector[id_vector %in% geneUniverse]
    
    validate(
      need(length(genesOfInterest) > 0, "No matching Gene IDs found in the gene universe.")
    )
    
    # Create a factor indicating genes of interest witBSU00240, BSU00260, BSU00280,BSU00290,BSU00300,BSU00310,BSU00320,BSU00330,BSU00340, BSU00350,BSU00360,BSU00370,BSU00380, BSU00390hin the gene universe
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    
    # Function to perform GO enrichment analysis
    performGOAnalysis <- function(ontology, gene2GO) {
      GOdata <- new("topGOdata",
                    description = "GO Enrichment Analysis",
                    ontology = ontology,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
      
      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
      
      # Retrieve and format the results
      allResults <- GenTable(GOdata,
                             classicFisher = resultFisher,
                             orderBy = "classicFisher",
                             ranksOf = "classicFisher",
                             topNodes = 10)
      
      # Adjust p-values
      adjP <- p.adjust(allResults$classicFisher, method = "fdr")
      allResults$FDR <- adjP
      
      return(allResults)
    }
    
    # Perform analysis for each ontology
    allResultsBP <- performGOAnalysis("BP", geneID2GObp)
    allResultsMF <- performGOAnalysis("MF", geneID2GOmf)
    allResultsCC <- performGOAnalysis("CC", geneID2GOcc)
    
    # Return a list of results
    list(BP = allResultsBP, MF = allResultsMF, CC = allResultsCC)
  })
  
  # Render the results in DataTables
  output$result_tableBP <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$BP, options = list(pageLength = 10))
  })
  
  output$result_tableMF <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$MF, options = list(pageLength = 10))
  })
  
  output$result_tableCC <- renderDT({
    req(analysis_results())
    datatable(analysis_results()$CC, options = list(pageLength = 10))
  })
  
  # Provide download functionality
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("GO_analysis_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(analysis_results())
      # Combine all results into one data frame
      all_results <- rbind(
        cbind(analysis_results()$BP, Ontology = "BP"),
        cbind(analysis_results()$MF, Ontology = "MF"),
        cbind(analysis_results()$CC, Ontology = "CC")
      )
      write.csv(all_results, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
