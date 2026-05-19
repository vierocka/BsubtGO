library(shiny)
library(DT)
library(topGO)

allGOwithConv <- read.table(
  "PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Build a named list mapping gene IDs (id_col) to GO terms (go_col).
# Rows with no GO annotations are silently skipped.
prepareGene2GO <- function(data, id_col, go_col) {
  universe <- as.character(data[, id_col])

  rows <- apply(data[, c(id_col, go_col)], 1, function(x) {
    terms <- unlist(strsplit(as.character(x[2]), split = " "))
    terms <- terms[nzchar(terms)]
    if (length(terms) == 0) return(NULL)
    list(id = x[1], terms = terms)
  })
  rows <- rows[!sapply(rows, is.null)]

  ids   <- sapply(rows, `[[`, "id")
  terms <- lapply(rows, `[[`, "terms")
  mapping <- setNames(terms, ids)
  mapping[names(mapping) %in% universe]
}

# Pre-compute mappings for both strains so startup is fast per-user
geneUniverse_bsub <- as.character(allGOwithConv$IDbsub)
geneUniverse_pg10 <- as.character(allGOwithConv$IDpg10)

gene2GO_bsub <- list(
  BP = prepareGene2GO(allGOwithConv, 2, 3),
  MF = prepareGene2GO(allGOwithConv, 2, 4),
  CC = prepareGene2GO(allGOwithConv, 2, 5)
)
gene2GO_pg10 <- list(
  BP = prepareGene2GO(allGOwithConv, 1, 3),
  MF = prepareGene2GO(allGOwithConv, 1, 4),
  CC = prepareGene2GO(allGOwithConv, 1, 5)
)

# topGO sometimes returns p-values as strings like "< 1e-30"
parseTopGOPval <- function(x) {
  suppressWarnings(as.numeric(sub("^<\\s*", "", x)))
}

ui <- fluidPage(
  titlePanel("Gene Ontology (GO) Term Analysis of Bacillus subtilis"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "strain",
        label = "Select strain of Bacillus subtilis:",
        choices = c("BSUB168", "PG10"),
        selected = "BSUB168"
      ),
      textAreaInput(
        inputId = "id_list",
        label = "Enter Gene IDs (comma- or newline-separated):",
        placeholder = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300",
        rows = 6
      ),
      actionButton("analyze", "Analyze", class = "btn-primary"),
      br(), br(),
      downloadButton("downloadData", "Download Results")
    ),
    mainPanel(
      h4("Biological Processes:"),
      br(),
      DTOutput("result_tableBP"),
      br(), br(),
      h4("Molecular Function:"),
      br(),
      DTOutput("result_tableMF"),
      br(), br(),
      h4("Cellular Components:"),
      br(),
      DTOutput("result_tableCC")
    )
  )
)

server <- function(input, output, session) {

  # Update placeholder text to match the selected strain
  observeEvent(input$strain, {
    if (input$strain == "PG10") {
      updateTextAreaInput(session, "id_list",
                          # PG10 always any
        placeholder = "ANY33920.1, ANY33921.1, ANY33922.1, ANY33923.1")
    } else {
      updateTextAreaInput(session, "id_list",
        placeholder = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300")
    }
  })

  analysis_results <- eventReactive(input$analyze, {
    req(input$id_list)

    id_vector <- unlist(strsplit(input$id_list, split = "[,\r\n]+"))
    id_vector <- trimws(id_vector)
    id_vector <- id_vector[nzchar(id_vector)]

    validate(need(length(id_vector) > 0, "Please enter at least one gene ID."))

    if (input$strain == "PG10") {
      universe <- geneUniverse_pg10
      g2go     <- gene2GO_pg10
    } else {
      universe <- geneUniverse_bsub
      g2go     <- gene2GO_bsub
    }

    genesOfInterest <- id_vector[id_vector %in% universe]
    validate(need(
      length(genesOfInterest) > 0,
      paste0(
        "No matching Gene IDs found for strain ", input$strain, ". ",
        "Make sure your IDs match the selected strain ",
        "(e.g. BSU00240 for BSUB168, ANY33920.1 for PG10)."
      )
    ))

    geneList <- factor(as.integer(universe %in% genesOfInterest))
    names(geneList) <- universe

    runOntology <- function(ontology) {
      GOdata <- new("topGOdata",
        description = "GO Enrichment Analysis",
        ontology    = ontology,
        allGenes    = geneList,
        annot       = annFUN.gene2GO,
        gene2GO     = g2go[[ontology]])

      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

      res <- GenTable(GOdata,
        classicFisher = resultFisher,
        orderBy       = "classicFisher",
        ranksOf       = "classicFisher",
        topNodes      = 20)

      res$FDR <- p.adjust(parseTopGOPval(res$classicFisher), method = "fdr")
      res
    }

    list(
      BP = runOntology("BP"),
      MF = runOntology("MF"),
      CC = runOntology("CC")
    )
  })

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

  output$downloadData <- downloadHandler(
    filename = function() paste0("GO_analysis_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(analysis_results())
      all_results <- rbind(
        cbind(analysis_results()$BP, Ontology = "BP"),
        cbind(analysis_results()$MF, Ontology = "MF"),
        cbind(analysis_results()$CC, Ontology = "CC")
      )
      write.csv(all_results, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
