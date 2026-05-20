# sensitivity_analysis.R
#
# Demonstrates how removing or swapping a single gene near the FDR 0.05
# boundary can dramatically reshape GO Biological Process enrichment results.
#
# Four scenarios (BSU168 IDs, BSU168 background, n = 4,185 genes):
#   A  full set     BSU00420–00500  (9 unique IDs; duplicate BSU00460 dropped)
#   B  drop one     A minus BSU00430
#   C  drop one     A minus BSU00500
#   D  swap one     A minus BSU00430, plus BSU00510
#
# Per set: top 10 BP terms by FDR (ascending).
# Union of all returned terms becomes the shared row axis of both panels.
#
# Figure — two panels, shared row labels:
#   Panel 1  presence / absence  — filled if term is in set's top 10
#   Panel 2  FDR value           — colour-coded FDR; grey where term absent
#
# Usage
#   cd BsubtGO/
#   Rscript sensitivity_analysis.R

library(ggplot2)
library(patchwork)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
bsub168_go <- read.table(
  "bsub168_go.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""
)

goTermNames <- local({
  df <- read.table("go_terms.tsv", sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, quote = "")
  setNames(df$Term, df$GO.ID)
})

# ---------------------------------------------------------------------------
# GO → gene index and Fisher test
# ---------------------------------------------------------------------------
makeGO2Gene <- function(data, id_col, go_col) {
  gene_ids   <- as.character(data[[id_col]])
  go_strings <- as.character(data[[go_col]])
  term_lists <- strsplit(go_strings, " ", fixed = TRUE)
  term_lists <- lapply(term_lists, function(x) x[nzchar(x)])
  split(rep(gene_ids, lengths(term_lists)), unlist(term_lists))
}

go2gene_bsub <- makeGO2Gene(bsub168_go, "IDbsub", "BiolProc")
universe_bsub <- as.character(bsub168_go$IDbsub)

top_BP <- function(genes, top_n = 10) {
  genes <- genes[genes %in% universe_bsub]
  n_u <- length(universe_bsub)
  n_s <- length(genes)
  if (n_s == 0L) return(NULL)

  rows <- lapply(names(go2gene_bsub), function(go_id) {
    ann   <- go2gene_bsub[[go_id]]
    n_sig <- sum(ann %in% genes)
    if (n_sig == 0L) return(NULL)
    n_ann <- length(ann)
    p <- fisher.test(
      matrix(c(n_sig, n_s - n_sig, n_ann - n_sig,
               n_u - n_s - n_ann + n_sig), nrow = 2),
      alternative = "greater"
    )$p.value
    data.frame(GO.ID = go_id, pvalue = p, stringsAsFactors = FALSE)
  })

  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0L) return(NULL)

  res      <- do.call(rbind, rows)
  res      <- res[order(res$pvalue), ]
  res$FDR  <- p.adjust(res$pvalue, method = "fdr")
  res$Term <- goTermNames[res$GO.ID]
  head(res, top_n)
}

# ---------------------------------------------------------------------------
# Four gene sets
# ---------------------------------------------------------------------------
sets <- list(
  "A\nfull\n(n=9)"         = c("BSU00420","BSU00430","BSU00440","BSU00450",
                                "BSU00460","BSU00470","BSU00480","BSU00490",
                                "BSU00500"),
  "B\n−BSU00430\n(n=8)"    = c("BSU00420","BSU00440","BSU00450","BSU00460",
                                "BSU00470","BSU00480","BSU00490","BSU00500"),
  "C\n−BSU00500\n(n=8)"   = c("BSU00420","BSU00430","BSU00440","BSU00450",
                                "BSU00460","BSU00470","BSU00480","BSU00490"),
  "D\n−BSU00430\n+BSU00510\n(n=9)" = c("BSU00420","BSU00440","BSU00450",
                                        "BSU00460","BSU00470","BSU00480",
                                        "BSU00490","BSU00500","BSU00510")
)

cat("Running GO BP enrichment for 4 gene sets...\n")
results <- lapply(sets, top_BP, top_n = 10)

for (nm in names(sets)) {
  r <- results[[nm]]
  cat(sprintf("  %-35s : %d terms returned\n",
              gsub("\n", " ", nm), if (is.null(r)) 0L else nrow(r)))
}

# ---------------------------------------------------------------------------
# Union of all GO terms; build tidy long-format data frame
# ---------------------------------------------------------------------------
all_ids <- unique(unlist(lapply(results, function(r) if (!is.null(r)) r$GO.ID)))

# Term labels: name (GO:ID), truncated to 45 chars
term_label <- function(go_id) {
  nm <- goTermNames[[go_id]]
  if (is.na(nm) || !nzchar(nm)) nm <- go_id
  lbl <- sprintf("%s (%s)", nm, go_id)
  if (nchar(lbl) > 50) lbl <- paste0(substr(lbl, 1, 47), "...")
  lbl
}
term_labels <- setNames(sapply(all_ids, term_label), all_ids)

# Order rows: most-frequent terms first, then by median FDR across sets
term_freq <- sapply(all_ids, function(id)
  sum(sapply(results, function(r) !is.null(r) && id %in% r$GO.ID)))
term_mfdr <- sapply(all_ids, function(id) {
  fdrs <- sapply(results, function(r) {
    if (is.null(r)) return(NA_real_)
    idx <- match(id, r$GO.ID)
    if (is.na(idx)) NA_real_ else r$FDR[idx]
  })
  median(fdrs, na.rm = TRUE)
})
row_order <- order(-term_freq, term_mfdr)
ordered_ids    <- all_ids[row_order]
ordered_labels <- term_labels[ordered_ids]

set_names_ordered <- names(sets)

df_long <- do.call(rbind, lapply(set_names_ordered, function(nm) {
  r <- results[[nm]]
  do.call(rbind, lapply(ordered_ids, function(id) {
    if (!is.null(r) && id %in% r$GO.ID) {
      fdr <- r$FDR[r$GO.ID == id]
      data.frame(term    = ordered_labels[id],
                 set     = nm,
                 present = TRUE,
                 fdr     = fdr,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(term    = ordered_labels[id],
                 set     = nm,
                 present = FALSE,
                 fdr     = NA_real_,
                 stringsAsFactors = FALSE)
    }
  }))
}))

df_long$term     <- factor(df_long$term, levels = rev(ordered_labels))
df_long$set      <- factor(df_long$set,  levels = set_names_ordered)
df_long$log10fdr <- -log10(df_long$fdr)   # NA stays NA; higher = more significant

# ---------------------------------------------------------------------------
# Panel 1 — presence / absence
# ---------------------------------------------------------------------------
p1 <- ggplot(df_long, aes(x = set, y = term, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "grey93"),
                    labels = c("TRUE" = "in top 10", "FALSE" = "absent"),
                    name   = NULL) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Presence / absence") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(size = 8.5, lineheight = 0.85),
        axis.text.y  = element_text(size = 7.5),
        panel.grid   = element_blank(),
        legend.position = "bottom",
        plot.title   = element_text(size = 10, face = "bold"))

# Short labels for panel 2 (full descriptions already in panel 1)
short_labels <- setNames(c("A", "B", "C", "D"), set_names_ordered)

# ---------------------------------------------------------------------------
# Panel 2 — FDR value
# ---------------------------------------------------------------------------
p2 <- ggplot(df_long, aes(x = set, y = term, fill = log10fdr)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low      = "#deebf7",
    high     = "#08306b",
    na.value = "grey93",
    name     = expression(-log[10](FDR)),
    guide    = guide_colorbar(barwidth = 6, barheight = 0.6)
  ) +
  scale_x_discrete(position = "top", labels = short_labels) +
  labs(x = NULL, y = NULL,
       title = "FDR value (dark = significant)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(size = 10, face = "bold"),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        legend.position = "bottom",
        plot.title   = element_text(size = 10, face = "bold"))

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
combined <- p1 + p2 +
  plot_layout(widths = c(1.6, 1)) +
  plot_annotation(
    title    = "Sensitivity of GO BP enrichment to single-gene changes",
    subtitle = sprintf(
      "BSU168 background (n = %d)  ·  top 10 terms per set  ·  %d unique terms shown\nSet A: %s",
      length(universe_bsub), length(all_ids),
      paste(sets[[1]], collapse = ", ")
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "grey40")
    )
  )

n_rows  <- length(all_ids)
fig_h   <- max(4, 1.2 + n_rows * 0.32)
ggsave("sensitivity_heatmap.png", combined,
       width = 9, height = fig_h, dpi = 150, limitsize = FALSE)
cat(sprintf("\nSaved: sensitivity_heatmap.png  (%d terms × 4 sets)\n",
            length(all_ids)))

# ---------------------------------------------------------------------------
# Save table
# ---------------------------------------------------------------------------
df_out <- do.call(rbind, lapply(names(sets), function(nm) {
  r <- results[[nm]]
  if (is.null(r)) return(NULL)
  r$set <- gsub("\n", " ", nm)
  r[, c("set", "GO.ID", "Term", "pvalue", "FDR")]
}))
write.csv(df_out, "sensitivity_analysis_results.csv", row.names = FALSE,
          quote = FALSE)
cat("Saved: sensitivity_analysis_results.csv\n")
