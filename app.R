library(topGO)
library(shiny)
library(DT)
library(clusterProfiler)

# setwd("~/Dropbox/Melih_AGMeier") # modify the path to the files

ui <- fluidPage(
  # Application title
  titlePanel("Gene Ontology (GO) Term Analysis of Bucillus subtilis"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        inputId = "id_list",
        label = "Enter Gene IDs separated by comma",
        placeholder = "BSU00240, BSU00260, BSU00280,BSU00290,BSU00300,BSU00310,BSU00320,BSU00330,BSU00340, BSU00350,BSU00360,BSU00370,BSU00380, BSU00390",
        rows = 6
      )
      ,
      selectInput(
        inputId = "keytype",
        label = "Select strain of Bacillus subtilis:",
        choices = c("PG10", "BSUB168"),
        selected = "PG10"
      ),
      #selectInput(
      #  inputId = "ontology",
      #  label = "Select Ontology:",
       # choices = c("Biological Process (BP)" = "BP",
      #              "Molecular Function (MF)" = "MF",
       #             "Cellular Component (CC)" = "CC"),
       # selected = "BP"
      #),
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
      h4("Cellular compartments:"),
      br(),
      DTOutput("result_tableCC"),
      
    )
  )
)


# Define server logic
server <- function(input, output){
  allGOwithConv <- read.table("PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv", sep="\t", header = TRUE)
  
  #### GENE LIST
  # NCBI: https://www.ncbi.nlm.nih.gov/nuccore/CP016788
  # download the gff3 file, protein fasta file
  # NCBI blastp script to compare Bsubtilis 168 and PG10 - only best hits
  # add GO terms from string-db to the equivalent Bsub168 ID: PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv
  # Read your gene universe (all PG10 genes)
  geneUniverse <- allGOwithConv$IDbsub
  # Read your genes of interest
  # example: ,"BSU00280","BSU00290","BSU00300","BSU00310","BSU00320","BSU00330","BSU00340","BSU00350","BSU00360","BSU00370","BSU00380","BSU00390","BSU00400"
  # genesOfInterest <- c("BSU00190","BSU00200","BSU00210","BSU00220","BSU00230","BSU00240","BSU00250","BSU00260","BSU00270")
  genesOfInterest <- eventReactive(input$analyze, {
    req(input$id_list)  # Ensure input is available
    
    # Split the input IDs into a vector
    id_vector <- unlist(strsplit(input$id_list, split = "[\r\n]+"))
    id_vector <- trimws(id_vector)  # Remove leading/trailing whitespace
    id_vector <- id_vector[nzchar(id_vector)]  # Remove empty strings
    
    validate(
      need(length(id_vector) > 0, "Please enter at least one gene ID.")
    ) })
    
  # 2. Create a gene-to-GO mapping
  
  # Assuming you have a CSV file with two columns: 'GeneID' and 'GOterm'
  geneID2GObiolp <- allGOwithConv[,c(2,3)]
  geneID2GOmolf <- allGOwithConv[,c(2,4)]
  geneID2GOcellc <- allGOwithConv[,c(2,5)]
  
  # Create a list where each element is a gene ID and its associated GO terms
  myExtensionFun <- function(x){
    GOlist <- unlist(strsplit(x[2], split = " "))
    IDarray <- rep(x[1], length(GOlist))
    ID_GO <- cbind(IDarray, GOlist)
    return(ID_GO)
  }
  
  result_listBP <- apply(geneID2GObiolp,1,myExtensionFun)
  result_listMF <- apply(geneID2GOmolf,1,myExtensionFun)
  result_listCC <- apply(geneID2GOcellc,1,myExtensionFun)
  
  ID_GO_df_BP <- do.call(rbind, result_listBP)
  ID_GO_df_BP <- as.data.frame(ID_GO_df_BP)
  
  ID_GO_df_MF <- do.call(rbind, result_listMF)
  ID_GO_df_MF <- as.data.frame(ID_GO_df_MF)
  
  ID_GO_df_CC <- do.call(rbind, result_listCC)
  ID_GO_df_CC <- as.data.frame(ID_GO_df_CC)
  
  # Optional: Rename the columns for clarity
  colnames(ID_GO_df_BP) <- c("GeneID", "BP")
  colnames(ID_GO_df_MF) <- c("GeneID", "MF")
  colnames(ID_GO_df_CC) <- c("GeneID", "CC")
  
  geneID2GObp <- by(ID_GO_df_BP$BP, ID_GO_df_BP$GeneID, function(x) as.character(x))
  geneID2GOmf <- by(ID_GO_df_MF$MF, ID_GO_df_MF$GeneID, function(x) as.character(x))
  geneID2GOcc <- by(ID_GO_df_CC$CC, ID_GO_df_CC$GeneID, function(x) as.character(x))
  
  # Ensure that the names of geneID2GO match your gene universe
  geneID2GObp <- geneID2GObp[names(geneID2GObp) %in% geneUniverse]
  geneID2GOmf <- geneID2GOmf[names(geneID2GOmf) %in% geneUniverse]
  geneID2GOcc <- geneID2GOcc[names(geneID2GOcc) %in% geneUniverse]
  
  # 3. Prepare the gene list for topGO
  
  # Create a factor indicating genes of interest within the gene universe
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # 4. Create a topGOdata object
  
  GOdatabp <- new("topGOdata",
                  description = "GO Enrichment Analysis",
                  ontology = "BP",  # You can also use "MF" or "CC"
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GObp)
  
  GOdatamf <- new("topGOdata",
                  description = "GO Enrichment Analysis",
                  ontology = "MF",  # You can also use "BP" or "CC"
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOmf)
  
  GOdatacc <- new("topGOdata",
                  description = "GO Enrichment Analysis",
                  ontology = "CC",  # You can also use "BP" or "MF"
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GOcc)
  
  # 5. Run the enrichment analysis
  
  resultFisherBP <- runTest(GOdatabp, algorithm = "classic", statistic = "fisher")
  resultFisherMF <- runTest(GOdatamf, algorithm = "classic", statistic = "fisher")
  resultFisherCC <- runTest(GOdatacc, algorithm = "classic", statistic = "fisher")
  
  # 6. Retrieve and display the results
  
  # You can adjust 'topNodes' to display more or fewer results
  allResultsBP <- GenTable(GOdatabp,
                           classicFisher = resultFisherBP,
                           orderBy = "classicFisher",
                           ranksOf = "classicFisher",
                           topNodes = 10)
  
  
  allResultsMF <- GenTable(GOdatamf,
                           classicFisher = resultFisherMF,
                           orderBy = "classicFisher",
                           ranksOf = "classicFisher",
                           topNodes = 10)
  
  
  allResultsCC <- GenTable(GOdatacc,
                           classicFisher = resultFisherCC,
                           orderBy = "classicFisher",
                           ranksOf = "classicFisher",
                           topNodes = 10)
  
  
  adjBP <- p.adjust(allResultsBP$classicFisher, method = "fdr")
  adjMF <- p.adjust(allResultsMF$classicFisher, method = "fdr")
  adjCC <- p.adjust(allResultsCC$classicFisher, method = "fdr")
  allResultsBP$FDR <- adjBP
  allResultsMF$FDR <- adjMF
  allResultsCC$FDR <- adjCC

  # Render the results in a DataTable
  output$result_tableBP <- renderDT({
    allResultsBP
  })
  
  output$result_tableMF <- renderDT({
    allResultsMF
  })
  
  output$result_tableCC <- renderDT({
    allResultsCC
  })
  # Provide download functionality
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("GO_analysis_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysis_results(), file, row.names = FALSE)
    }
  )

  }

# Run the application
shinyApp(ui = ui, server = server)

