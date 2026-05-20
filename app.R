library(shiny)
library(DT)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

# PG10 <-> BSU168 mapping table (2,683 rows; one row per PG10 gene with a
# high-confidence BLASTP best-hit BSU168 homolog, e-value < 1e-6).
# GO terms are BSU168-derived annotations transferred via BLASTP.
mapping <- read.table(
  "PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Full BSU168 GO table: all 4,185 proteins, GO terms from STRING-db v12.0.
bsub168_go <- read.table(
  "bsub168_go.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""
)

# GO term ID -> human-readable name
goTermNames <- local({
  df <- read.table("go_terms.tsv", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "")
  setNames(df$Term, df$GO.ID)
})

# ---------------------------------------------------------------------------
# Build inverted GO->gene indices (vectorised; runs once at startup)
# ---------------------------------------------------------------------------

prepareGO2Gene <- function(data, id_col, go_col) {
  gene_ids   <- as.character(data[[id_col]])
  go_strings <- as.character(data[[go_col]])
  term_lists <- strsplit(go_strings, " ", fixed = TRUE)
  term_lists <- lapply(term_lists, function(x) x[nzchar(x)])
  genes_expanded <- rep(gene_ids, lengths(term_lists))
  terms_expanded <- unlist(term_lists)
  split(genes_expanded, terms_expanded)
}

# Mode 1 – PG10 IDs, PG10 background
go2gene_pg10 <- list(
  BP = prepareGO2Gene(mapping, "IDpg10", "BiolProc"),
  MF = prepareGO2Gene(mapping, "IDpg10", "MolFun"),
  CC = prepareGO2Gene(mapping, "IDpg10", "CellComp")
)
universe_pg10 <- as.character(mapping$IDpg10)   # 2,683 PG10 gene IDs

# Mode 2 – BSU168 IDs entered, but PG10 background (same indices as Mode 1;
# user IDs are converted BSU->PG10 before the test)
# (universe_pg10 and go2gene_pg10 are reused)

# Mode 3 – BSU168 IDs, full BSU168 background (4,185 genes)
go2gene_bsub <- list(
  BP = prepareGO2Gene(bsub168_go, "IDbsub", "BiolProc"),
  MF = prepareGO2Gene(bsub168_go, "IDbsub", "MolFun"),
  CC = prepareGO2Gene(bsub168_go, "IDbsub", "CellComp")
)
universe_bsub <- as.character(bsub168_go$IDbsub)   # 4,185 BSU168 gene IDs

# ---------------------------------------------------------------------------
# Fisher exact test for GO over-representation
# Equivalent to the classic/Fisher algorithm in topGO.
# ---------------------------------------------------------------------------

runGOFisher <- function(go2gene, genes_of_interest, universe) {
  n_universe <- length(universe)
  n_selected <- length(genes_of_interest)

  rows <- lapply(names(go2gene), function(go_id) {
    annotated   <- go2gene[[go_id]]
    n_annotated <- length(annotated)
    n_sig       <- sum(annotated %in% genes_of_interest)
    if (n_sig == 0) return(NULL)

    mat <- matrix(c(
      n_sig,
      n_selected  - n_sig,
      n_annotated - n_sig,
      n_universe  - n_selected - (n_annotated - n_sig)
    ), nrow = 2)

    p <- fisher.test(mat, alternative = "greater")$p.value

    data.frame(
      GO.ID       = go_id,
      Term        = unname(goTermNames[go_id]),
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

  res     <- do.call(rbind, rows)
  res     <- res[order(res$pvalue), ]
  res$FDR <- signif(p.adjust(res$pvalue, method = "fdr"), 4)
  head(res, 20)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

mode_choices <- c(
  "PG10 genes  |  PG10 background (2,683 genes)"           = "pg10",
  "PG10 genes via BSU168 IDs  |  PG10 background (2,683 genes)" = "pg10_bsu",
  "BSU168 genes  |  full BSU168 background (4,185 genes)"  = "bsub168"
)

about_tab <- tabPanel(
  "About & Methods",
  fluidRow(column(10, offset = 1,
    h3("BsubtilisGO — background"),
    p("GO enrichment analysis asks whether any biological functions are
      statistically over-represented in a gene list compared to what would be
      expected by chance. The answer depends critically on the ",
      strong("background gene set"), " (the universe): every gene that
      could, in principle, have appeared in the list."),
    p("Most online tools (g:Profiler, DAVID, STRING web interface) only carry
      ", em("B. subtilis"), " strain ", strong("168"), " in their databases.
      Analysing an experiment done in ", strong("PG10"), " — a reduced-genome
      strain with ~2,769 protein-coding genes — with a 168 background
      introduces two errors:"),
    tags$ul(
      tags$li("Genes present in 168 but absent from PG10 dilute enrichment
               signals and cause false negatives."),
      tags$li("PG10-specific genes have no 168 equivalent and are silently
               dropped, distorting the statistics.")
    ),
    p("BsubtilisGO addresses this by building strain-matched backgrounds.
      GO annotations are transferred from the well-annotated strain 168 to
      PG10 via BLASTP best-hit homology, and the enrichment test is run
      against only the genes that actually exist in the strain of interest."),

    hr(),
    h3("The three analysis modes"),

    h4("1. PG10 genes | PG10 background"),
    p("Input: PG10 protein IDs (", code("ANY*.1"), " format, e.g. ",
      code("ANY33920.1"), ")."),
    p("Background: all 2,683 PG10 proteins for which a high-confidence BSU168
      BLASTP best-hit exists (e-value < 1×10⁻⁶). GO annotations are inherited
      from the matched BSU168 homolog."),
    p("Use this mode when your gene list comes directly from a PG10 experiment
      and you want a fully PG10-native analysis."),

    h4("2. PG10 genes via BSU168 IDs | PG10 background"),
    p("Input: BSU168 locus tags (", code("BSU*"), " format, e.g. ",
      code("BSU00240"), ") for genes you know are present in PG10."),
    p("The app converts each BSU168 ID to its PG10 equivalent(s) via the
      BLASTP mapping table, then runs the exact same test as mode 1 with the
      PG10 background."),
    p("Use this mode when your gene list is expressed in BSU168 locus tags
      (e.g. from literature) but the experiment was done in PG10."),

    h4("3. BSU168 genes | full BSU168 background"),
    p("Input: BSU168 locus tags (", code("BSU*"), " format)."),
    p("Background: the full BSU168 proteome of 4,185 proteins with GO
      annotations loaded directly from STRING-db v12.0 — not filtered through
      the PG10 mapping. Useful for demonstrating how a mismatched background
      inflates or deflates enrichment compared to mode 1."),
    p("Use this mode for experiments done in strain 168, or when you want the
      complete BSU168 reference without any PG10-driven truncation."),

    hr(),
    h3("Statistical method"),
    p("Each GO term is tested with a one-sided Fisher exact test
      (over-representation, ", code('alternative = "greater"'), "). The
      2 × 2 contingency table counts genes annotated / not annotated to the
      term, within / outside the submitted list. This is equivalent to the ",
      strong("classic/Fisher algorithm"), " in the topGO package ",
      "(Alexa A & Rahnenführer J, 2023, Bioconductor)."),
    p("Multiple-testing correction: Benjamini–Hochberg FDR, applied to the
      top 20 results per ontology. Because FDR is calculated on the truncated
      list rather than all tested terms, treat borderline values with caution."),

    hr(),
    h3("Data sources"),
    tags$ul(
      tags$li(em("B. subtilis"), " 168 proteome and GO annotations: STRING-db
               v12.0, species 224308."),
      tags$li("PG10 proteome: NCBI accession CP016788."),
      tags$li("PG10–BSU168 homology: BLASTP best-hit, e-value < 1×10⁻⁶
               (86 low-confidence hits from the raw output were removed).")
    )
  ))
)

ui <- fluidPage(
  titlePanel("Gene Ontology (GO) Enrichment Analysis — Bacillus subtilis"),
  tabsetPanel(
    tabPanel(
      "Analysis",
      br(),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId  = "mode",
            label    = "Analysis mode:",
            choices  = mode_choices,
            selected = "pg10"
          ),
          uiOutput("mode_info"),
          br(),
          uiOutput("id_input_ui"),
          actionButton("analyze", "Analyze", class = "btn-primary"),
          br(), br(),
          downloadButton("downloadData", "Download Results")
        ),
        mainPanel(
          h4("Biological Process:"),
          DTOutput("result_tableBP"),
          br(),
          h4("Molecular Function:"),
          DTOutput("result_tableMF"),
          br(),
          h4("Cellular Component:"),
          DTOutput("result_tableCC")
        )
      )
    ),
    about_tab
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # Dynamic help text under the mode selector
  output$mode_info <- renderUI({
    msg <- switch(input$mode,
      pg10 = tags$small(style = "color:#555;",
        "Enter PG10 protein IDs (", code("ANY*.1"), "). ",
        "Background: 2,683 PG10 genes with a high-confidence BSU168 homolog ",
        "(BLASTP best-hit, e-value < 1×10⁻⁶)."
      ),
      pg10_bsu = tags$small(style = "color:#555;",
        "Enter BSU168 locus tags (", code("BSU*"), ") for genes present in PG10. ",
        "The app maps them to PG10 IDs and tests against the PG10 background ",
        "(2,683 genes)."
      ),
      bsub168 = tags$small(style = "color:#555;",
        "Enter BSU168 locus tags (", code("BSU*"), "). ",
        "Background: full BSU168 proteome, 4,185 genes, GO from STRING-db v12.0."
      )
    )
    msg
  })

  # Dynamic text area with mode-appropriate placeholder and label
  output$id_input_ui <- renderUI({
    cfg <- switch(input$mode,
      pg10     = list(ph = "ANY33920.1, ANY33921.1, ANY33922.1, ANY33923.1",
                      lb = "PG10 protein IDs (comma- or newline-separated):"),
      pg10_bsu = list(ph = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300",
                      lb = "BSU168 locus tags for PG10 genes (comma- or newline-separated):"),
      bsub168  = list(ph = "BSU00240, BSU00260, BSU00280, BSU00290, BSU00300",
                      lb = "BSU168 locus tags (comma- or newline-separated):")
    )
    textAreaInput("id_list", label = cfg$lb, placeholder = cfg$ph, rows = 6)
  })

  analysis_results <- eventReactive(input$analyze, {
    req(input$id_list)

    raw <- unlist(strsplit(input$id_list, "[,\r\n]+"))
    raw <- trimws(raw)
    raw <- raw[nzchar(raw)]
    validate(need(length(raw) > 0, "Please enter at least one gene ID."))

    if (input$mode == "pg10") {
      universe <- universe_pg10
      g2go     <- go2gene_pg10
      genes    <- raw[raw %in% universe]
      validate(need(length(genes) > 0,
        paste0("No PG10 IDs recognised. Use ANY*.1 format, e.g. ANY33920.1.")))

    } else if (input$mode == "pg10_bsu") {
      # Convert BSU168 IDs to PG10 IDs via mapping table
      matched  <- mapping[mapping$IDbsub %in% raw, ]
      validate(need(nrow(matched) > 0,
        paste0("None of the entered BSU168 IDs were found in the PG10 mapping. ",
               "Check that your IDs use BSU* format and are present in PG10.")))
      genes    <- unique(as.character(matched$IDpg10))
      universe <- universe_pg10
      g2go     <- go2gene_pg10
      # Report how many BSU IDs were converted
      n_entered   <- length(unique(raw))
      n_converted <- length(unique(matched$IDbsub))
      if (n_converted < n_entered) {
        showNotification(
          paste0(n_entered - n_converted, " of ", n_entered,
                 " BSU168 IDs had no PG10 equivalent and were skipped."),
          type = "warning", duration = 8
        )
      }

    } else {
      universe <- universe_bsub
      g2go     <- go2gene_bsub
      genes    <- raw[raw %in% universe]
      validate(need(length(genes) > 0,
        paste0("No BSU168 IDs recognised. Use BSU* format, e.g. BSU00240.")))
    }

    list(
      BP    = runGOFisher(g2go$BP, genes, universe),
      MF    = runGOFisher(g2go$MF, genes, universe),
      CC    = runGOFisher(g2go$CC, genes, universe),
      mode  = input$mode,
      n_in  = length(raw),
      n_matched = length(genes)
    )
  })

  render_table <- function(ont) {
    renderDT({
      req(analysis_results())
      res <- analysis_results()[[ont]]
      datatable(res, options = list(pageLength = 10),
                caption = if (nrow(res) == 0) "No significant terms found." else NULL)
    })
  }
  output$result_tableBP <- render_table("BP")
  output$result_tableMF <- render_table("MF")
  output$result_tableCC <- render_table("CC")

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("GO_results_", input$mode, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(analysis_results())
      r <- analysis_results()
      all_res <- rbind(
        cbind(r$BP, Ontology = "BP"),
        cbind(r$MF, Ontology = "MF"),
        cbind(r$CC, Ontology = "CC")
      )
      write.csv(all_res, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
