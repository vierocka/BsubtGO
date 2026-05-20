# background_comparison.R
#
# Demonstrates the effect of using a mismatched background (full BSU168,
# 4,185 genes) versus the correct strain-matched background (PG10, 2,683
# genes) on GO Biological Process enrichment FDR.
#
# Design
# ------
# 1. Sort the 2,678 unique BSU IDs present in the PG10 mapping into an array.
# 2. Sample 1,000 random starting positions (sliding windows of 10 consecutive
#    IDs). Consecutive BSU IDs cluster by operon/function, so some windows hit
#    real functional groups and produce low FDRs — unlike random 10-gene sets
#    which almost always return FDR ≈ 1.
# 3. For each window run two one-sided Fisher exact tests (Biological Process
#    only, equivalent to the classic/Fisher algorithm in topGO):
#      Correct : convert BSU IDs → IDpg10 via mapping table;
#                universe = 2,683 PG10 genes
#      Wrong   : use BSU IDs directly;
#                universe = 4,185 full BSU168 genes
# 4. Record the minimum FDR across all returned BP terms (= the top hit's FDR,
#    i.e. the LOWEST = most significant value). Windows with no annotated genes
#    return NA, recorded as 1.0.
# 5. Output: density plot + scatter plot comparing the two FDR distributions.
#
# Note on FDR behaviour
# ---------------------
# runGOFisher (app.R) has NO hard-coded FDR cutoff. Any term where at least
# one submitted gene is annotated (n_sig > 0) receives a Fisher p-value and a
# BH-adjusted FDR. Non-enriched windows return FDR ≈ 1; only windows where
# every gene is unannotated return NULL (treated as NA/1.0 here).
#
# Usage
# -----
#   cd BsubtGO/
#   Rscript background_comparison.R

library(ggplot2)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
mapping <- read.table(
  "PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

bsub168_go <- read.table(
  "bsub168_go.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = ""
)

# ---------------------------------------------------------------------------
# Build GO→gene inverted indices (Biological Process only)
# ---------------------------------------------------------------------------
makeGO2Gene <- function(data, id_col, go_col) {
  gene_ids   <- as.character(data[[id_col]])
  go_strings <- as.character(data[[go_col]])
  term_lists <- strsplit(go_strings, " ", fixed = TRUE)
  term_lists <- lapply(term_lists, function(x) x[nzchar(x)])
  split(rep(gene_ids, lengths(term_lists)), unlist(term_lists))
}

go2gene_pg10 <- makeGO2Gene(mapping,    "IDpg10", "BiolProc")
go2gene_bsub <- makeGO2Gene(bsub168_go, "IDbsub", "BiolProc")

universe_pg10 <- as.character(mapping$IDpg10)       # 2,683 genes
universe_bsub <- as.character(bsub168_go$IDbsub)    # 4,185 genes

cat(sprintf("PG10 universe : %d genes  |  %d BP GO terms\n",
            length(universe_pg10), length(go2gene_pg10)))
cat(sprintf("BSU168 universe: %d genes  |  %d BP GO terms\n",
            length(universe_bsub), length(go2gene_bsub)))

# ---------------------------------------------------------------------------
# minFDR_BP: run Fisher over all BP terms; return the minimum (best) FDR.
# Returns NA when no gene in the set is annotated to any BP term.
# ---------------------------------------------------------------------------
minFDR_BP <- function(go2gene, genes, universe) {
  n_universe <- length(universe)
  n_selected <- length(genes)
  if (n_selected == 0L) return(NA_real_)

  pvals <- vapply(go2gene, function(annotated) {
    n_sig <- sum(annotated %in% genes)
    if (n_sig == 0L) return(NA_real_)
    n_ann <- length(annotated)
    fisher.test(
      matrix(c(n_sig,
               n_selected - n_sig,
               n_ann      - n_sig,
               n_universe - n_selected - n_ann + n_sig),
             nrow = 2),
      alternative = "greater"
    )$p.value
  }, numeric(1L))

  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0L) return(NA_real_)
  min(p.adjust(pvals, method = "fdr"))
}

# ---------------------------------------------------------------------------
# Sliding-window sampling
# ---------------------------------------------------------------------------
bsu_sorted  <- sort(unique(mapping$IDbsub))   # 2,678 BSU IDs, ordered
window_size <- 10L
n_windows   <- length(bsu_sorted) - window_size + 1L   # 2,669

set.seed(42)
starts <- sample(n_windows, 1000L, replace = FALSE)

cat(sprintf("\nRunning 1,000 windows (size %d) × 2 backgrounds...\n", window_size))

fdr_pg10 <- numeric(length(starts))
fdr_bsub <- numeric(length(starts))

for (i in seq_along(starts)) {
  bsu_win <- bsu_sorted[starts[i]:(starts[i] + window_size - 1L)]

  # Correct background: PG10 IDpg10 IDs vs PG10 universe
  pg10_genes  <- mapping$IDpg10[mapping$IDbsub %in% bsu_win]
  fdr_pg10[i] <- minFDR_BP(go2gene_pg10, pg10_genes, universe_pg10)

  # Wrong background: BSU IDs directly vs full BSU168 universe
  bsu_genes   <- bsu_win[bsu_win %in% universe_bsub]
  fdr_bsub[i] <- minFDR_BP(go2gene_bsub, bsu_genes, universe_bsub)

  if (i %% 200L == 0L) cat(sprintf("  %d / %d\n", i, length(starts)))
}

# NA → 1.0 (no annotations = no enrichment possible)
fdr_pg10[is.na(fdr_pg10)] <- 1
fdr_bsub[is.na(fdr_bsub)] <- 1

cat(sprintf(
  "\nWindows with min-FDR < 0.05 :  PG10 = %3d  |  BSU168 = %3d\n",
  sum(fdr_pg10 < 0.05), sum(fdr_bsub < 0.05)
))
cat(sprintf(
  "Windows with min-FDR < 0.20 :  PG10 = %3d  |  BSU168 = %3d\n",
  sum(fdr_pg10 < 0.20), sum(fdr_bsub < 0.20)
))
cat(sprintf(
  "Windows with min-FDR = 1.00 :  PG10 = %3d  |  BSU168 = %3d\n",
  sum(fdr_pg10 == 1),   sum(fdr_bsub == 1)
))

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

# Per-window data
df_windows <- data.frame(
  window_start = bsu_sorted[starts],
  window_end   = bsu_sorted[starts + window_size - 1L],
  minFDR_PG10  = fdr_pg10,
  minFDR_BSU168 = fdr_bsub,
  significant  = ifelse(fdr_pg10 < 0.05 & fdr_bsub < 0.05, "Both",
                 ifelse(fdr_pg10 < 0.05 & fdr_bsub >= 0.05, "PG10_only",
                 ifelse(fdr_pg10 >= 0.05 & fdr_bsub < 0.05, "BSU168_only",
                                                              "Neither")))
)
write.csv(df_windows, "background_comparison_windows.csv", row.names = FALSE,
          quote = FALSE)

# Quadrant summary
df_summary <- data.frame(
  category    = c("Both significant (FDR < 0.05)",
                  "PG10 only — false negatives with BSU168",
                  "BSU168 only — false positives with BSU168",
                  "Neither significant"),
  n_windows   = c(sum(fdr_pg10 < 0.05 & fdr_bsub < 0.05),
                  sum(fdr_pg10 < 0.05 & fdr_bsub >= 0.05),
                  sum(fdr_pg10 >= 0.05 & fdr_bsub < 0.05),
                  sum(fdr_pg10 >= 0.05 & fdr_bsub >= 0.05)),
  pct         = round(c(sum(fdr_pg10 < 0.05 & fdr_bsub < 0.05),
                        sum(fdr_pg10 < 0.05 & fdr_bsub >= 0.05),
                        sum(fdr_pg10 >= 0.05 & fdr_bsub < 0.05),
                        sum(fdr_pg10 >= 0.05 & fdr_bsub >= 0.05)) / length(starts) * 100, 1)
)
write.csv(df_summary, "background_comparison_summary.csv", row.names = FALSE,
          quote = FALSE)

cat("\nSaved: background_comparison_windows.csv\n")
cat(  "       background_comparison_summary.csv\n")

# ---------------------------------------------------------------------------
# Plot 1: overlapping density — FDR distributions
# ---------------------------------------------------------------------------
col_pg10 <- "#2166ac"   # blue  = correct
col_bsub <- "#d73027"   # red   = wrong

df_dens <- rbind(
  data.frame(min_FDR = fdr_pg10, Background = "PG10 correct (n = 2,683)"),
  data.frame(min_FDR = fdr_bsub, Background = "BSU168 wrong  (n = 4,185)")
)
df_dens$Background <- factor(df_dens$Background,
  levels = c("PG10 correct (n = 2,683)", "BSU168 wrong  (n = 4,185)"))

p1 <- ggplot(df_dens, aes(x = min_FDR, fill = Background, colour = Background)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(values   = c("PG10 correct (n = 2,683)" = col_pg10,
                                  "BSU168 wrong  (n = 4,185)"  = col_bsub)) +
  scale_colour_manual(values = c("PG10 correct (n = 2,683)" = col_pg10,
                                  "BSU168 wrong  (n = 4,185)"  = col_bsub)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.07, y = Inf, label = "FDR 0.05",
           hjust = 0, vjust = 1.4, size = 3.2, colour = "grey40") +
  labs(
    title    = "Background set effect on GO Biological Process enrichment",
    subtitle = sprintf("1,000 sliding windows of %d consecutive PG10-mapped BSU168 genes", window_size),
    x        = "Minimum FDR per window (lowest = most significant BP term)",
    y        = "Density",
    fill = NULL, colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.68, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---------------------------------------------------------------------------
# Plot 2: scatter — per-window PG10 FDR vs BSU168 FDR
# ---------------------------------------------------------------------------
df_scat <- data.frame(
  pg10  = fdr_pg10,
  bsub  = fdr_bsub,
  # quadrant labels for colouring
  quad  = ifelse(fdr_pg10 < 0.05 & fdr_bsub >= 0.05, "PG10 only",
          ifelse(fdr_pg10 >= 0.05 & fdr_bsub < 0.05, "BSU168 only",
          ifelse(fdr_pg10 < 0.05 & fdr_bsub < 0.05,  "Both",
                                                       "Neither")))
)
df_scat$quad <- factor(df_scat$quad,
  levels = c("Both", "PG10 only", "BSU168 only", "Neither"))

quad_cols <- c("Both"        = "#4dac26",
               "PG10 only"   = col_pg10,
               "BSU168 only" = col_bsub,
               "Neither"     = "grey70")

p2 <- ggplot(df_scat, aes(x = pg10, y = bsub, colour = quad)) +
  geom_point(alpha = 0.5, size = 1.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
  geom_hline(yintercept = 0.05, linetype = "dotted", colour = col_bsub) +
  geom_vline(xintercept = 0.05, linetype = "dotted", colour = col_pg10) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  scale_colour_manual(values = quad_cols) +
  labs(
    title    = "Per-window FDR: correct vs wrong background",
    subtitle = "Dashed diagonal = perfect agreement; dotted lines = FDR 0.05",
    x        = "Min FDR — PG10 background (correct)",
    y        = "Min FDR — BSU168 background (wrong)",
    colour   = "Significant at FDR 0.05"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.75, 0.18),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
ggsave("background_effect_density.png", p1, width = 7, height = 4.5, dpi = 150)
ggsave("background_effect_scatter.png", p2, width = 5.5, height = 5.5, dpi = 150)
cat("\nSaved: background_effect_density.png\n")
cat(  "       background_effect_scatter.png\n")
