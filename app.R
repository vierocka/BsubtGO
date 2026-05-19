library(shiny)
library(DT)

allGOwithConv <- read.table(
  "PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

goTermNames <- local({
  df <- read.table("go_terms.tsv", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "")
  setNames(df$Term, df$GO.ID)
})

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

# Pre-compute mappings for both strains
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

# Fisher exact test for GO over-representation (one-sided, "greater").
# Equivalent to topGO classic/Fisher but with no Bioconductor dependency.
runGOFisher <- function(gene2go, genes_of_interest, universe) {
  n_universe <- length(universe)
  n_selected <- length(genes_of_interest)

  all_terms <- unique(unlist(gene2go))

  rows <- lapply(all_terms, function(go_id) {
    annotated     <- names(gene2go)[vapply(gene2go, function(x) go_id %in% x, logical(1))]
    n_annotated   <- length(annotated)
    n_sig         <- sum(genes_of_interest %in% annotated)
    if (n_sig == 0) return(NULL)

    mat <- matrix(c(
      n_sig,
      n_selected - n_sig,
      n_annotated - n_sig,
      n_universe - n_selected - (n_annotated - n_sig)
    ), nrow = 2)

    p <- fisher.test(mat, alternative = "greater")$p.value

    data.frame(
      GO.ID       = go_id,
      Term        = unname(goTermNames[go_id] %||% go_id),
      Annotated   = n_annotated,
      Significant = n_sig,
      Expected    = round(n_selected * n_annotated / n_universe, 2),
      pvalue      = signif(p, 4),
      stringsAsFactors = FALSE
    )
  })

  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(data.frame(
    GO.ID = character(), Term = character(), Annotated = integer(),
    Significant = integer(), Expected = numeric(),
    pvalue = numeric(), FDR = numeric()
  ))

  res      <- do.call(rbind, rows)
  res      <- res[order(res$pvalue), ]
  res$FDR  <- signif(p.adjust(res$pvalue, method = "fdr"), 4)
  head(res, 20)
}

`%||%` <- function(a, b) if (!is.na(a) && nzchar(a)) a else b

ui <- fluidPage(
  titlePanel("Gene Ontology (GO) Term Analysis of Bacillus subtilis"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "strain",
        label   = "Select strain of Bacillus subtilis:",
        choices  = c("BSUB168", "PG10"),
        selected = "BSUB168"
      ),
      textAreaInput(
        inputId     = "id_list",
        label       = "Enter Gene IDs (comma- or newline-separated):",
        placeholder = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300",
        rows        = 6
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

  observeEvent(input$strain, {
    if (input$strain == "PG10") {
      updateTextAreaInput(session, "id_list",
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

    list(
      BP = runGOFisher(g2go$BP, genesOfInterest, universe),
      MF = runGOFisher(g2go$MF, genesOfInterest, universe),
      CC = runGOFisher(g2go$CC, genesOfInterest, universe)
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
    content  = function(file) {
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
