import json, uuid, pathlib

def cell_id():
    return str(uuid.uuid4())[:8]

def md(text):
    return {"cell_type": "markdown", "id": cell_id(), "metadata": {}, "source": text}

def code(text):
    return {"cell_type": "code", "execution_count": None, "id": cell_id(),
            "metadata": {}, "outputs": [], "source": text}

cells = []

# ── TITLE ───────────────────────────────────────────────────────────────────
cells.append(md(
"""# GO Enrichment: Background Effects and Single-Gene Sensitivity
## *Bacillus subtilis* PG10 strain — Python re-implementation

**Based on:** `background_comparison.R` and `sensitivity_analysis.R`
(BsubtGO project, AG Meier)

---

This notebook is a step-by-step Python translation of two R analyses that
illustrate two fundamental pitfalls in GO enrichment testing:

1. **Background comparison** — using the wrong reference universe (e.g. the
   full *B. subtilis* BSU168 proteome instead of the actual PG10 strain
   subset) systematically inflates or deflates FDR estimates.

2. **Single-gene sensitivity** — removing or swapping just one gene near the
   FDR 0.05 boundary can dramatically change which GO Biological Process
   terms appear significant.

Both analyses use a **from-scratch Fisher exact test + Benjamini–Hochberg FDR**
implementation — no topGO, goatools, or other GO libraries required.
This makes the statistical machinery fully transparent.

**Data files** (all in the same folder as this notebook):

| File | Content |
|---|---|
| `PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv` | PG10↔BSU168 ID mapping + GO annotations |
| `bsub168_go.tsv` | Full BSU168 strain GO annotations (4,185 genes) |
| `go_terms.tsv` | GO term ID → human-readable name lookup |
"""
))

# ── IMPORTS ─────────────────────────────────────────────────────────────────
cells.append(md("## 1. Setup"))

cells.append(code(
"""import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from IPython.display import Image, display

plt.rcParams["figure.dpi"] = 120
plt.rcParams["font.size"]  = 11
print("Packages loaded.")
"""
))

# ── DATA LOADING ─────────────────────────────────────────────────────────────
cells.append(md(
"""## 2. Data Loading

### Two gene universes

This study compares two *B. subtilis* strains:

- **BSU168** — the standard laboratory reference strain; fully annotated;
  4,185 protein-coding genes with GO annotations.
- **PG10** — an engineering strain derived from BSU168 by large chromosomal
  deletions; only 2,683 genes remain.

The mapping file (`PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv`) links each
PG10 gene (`IDpg10`) to its BSU168 counterpart (`IDbsub`) and carries the
GO annotations split into three columns: `BiolProc`, `MolFun`, `CellComp`.
Each cell contains space-separated GO term IDs (e.g. `GO:0006260 GO:0051052`).

The third file (`go_terms.tsv`) maps GO IDs to human-readable names —
used only in the sensitivity analysis to label heatmap rows.
"""
))

cells.append(code(
"""# Load PG10 mapping (2,683 PG10 genes, each with its BSU168 counterpart + GO terms)
mapping = pd.read_csv("PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv", sep="\\t",
                      dtype=str).fillna("")

# Load full BSU168 GO annotations (4,185 genes)
bsub168 = pd.read_csv("bsub168_go.tsv", sep="\\t", dtype=str).fillna("")

# Load GO term name lookup
go_names = pd.read_csv("go_terms.tsv", sep="\\t", dtype=str)
go_names = go_names.set_index("GO.ID")["Term"].to_dict()

print(f"PG10 mapping rows  : {len(mapping):,}")
print(f"BSU168 rows        : {len(bsub168):,}")
print(f"GO term names loaded: {len(go_names):,}")
print()
print("Mapping columns:", mapping.columns.tolist())
print(mapping.head(2).to_string())
"""
))

# ── GO INDEX ─────────────────────────────────────────────────────────────────
cells.append(md(
"""## 3. Core Functions

### 3.1 Building an inverted GO→gene index

GO enrichment requires knowing, for every GO term, which genes are annotated
to it.  We build an **inverted index** — a dictionary
`{GO_term_ID: [gene_id_1, gene_id_2, ...]}` — once per universe, then reuse
it for every Fisher test.

```
Input  :  gene_id    GO terms (space-separated in one cell)
          geneA      GO:0006260 GO:0051052
          geneB      GO:0006260
          geneC      GO:0051052

Output :  { "GO:0006260": ["geneA", "geneB"],
            "GO:0051052": ["geneA", "geneC"] }
```
"""
))

cells.append(code(
"""def make_go2gene(df, id_col, go_col):
    \"\"\"
    Build inverted index GO term → list of gene IDs.

    Parameters
    ----------
    df     : pd.DataFrame with at least id_col and go_col
    id_col : column name for gene IDs
    go_col : column name containing space-separated GO term strings

    Returns
    -------
    dict  {go_term_id: [gene_id, ...]}
    \"\"\"
    go2gene = {}
    for _, row in df.iterrows():
        gene_id  = str(row[id_col])
        terms    = [t for t in str(row[go_col]).split() if t.startswith("GO:")]
        for term in terms:
            go2gene.setdefault(term, []).append(gene_id)
    return go2gene

# Build Biological Process indices for both universes
go2gene_pg10  = make_go2gene(mapping,  "IDpg10",  "BiolProc")
go2gene_bsub  = make_go2gene(bsub168,  "IDbsub",  "BiolProc")

universe_pg10 = mapping["IDpg10"].tolist()    # 2,683 PG10 genes
universe_bsub = bsub168["IDbsub"].tolist()    # 4,185 BSU168 genes

print(f"PG10  universe : {len(universe_pg10):,} genes | {len(go2gene_pg10):,} BP GO terms")
print(f"BSU168 universe: {len(universe_bsub):,} genes | {len(go2gene_bsub):,} BP GO terms")
"""
))

cells.append(md(
"""### 3.2 Fisher exact test for GO enrichment

For each GO term we fill a 2×2 contingency table and apply a one-sided
Fisher exact test (alternative = "greater" = overrepresentation):

```
                    In gene set     Not in gene set     Total
Annotated to term      n_sig         n_ann − n_sig      n_ann
Not annotated       n_sel − n_sig   N−n_sel−n_ann+n_sig  N−n_ann
Total               n_selected       N − n_selected        N
```

where:
- `N`          = universe size (total genes with GO annotations)
- `n_selected` = number of genes in the submitted set
- `n_ann`      = number of universe genes annotated to this GO term
- `n_sig`      = overlap: genes in set AND annotated to this term

After computing p-values across all testable terms (those with `n_sig ≥ 1`),
Benjamini–Hochberg FDR correction is applied.  We return the **minimum FDR**
across all terms — i.e. the best (lowest) p-value after correction, which is
what gets plotted in the background comparison analysis.
"""
))

cells.append(code(
"""def min_fdr_bp(go2gene, genes, universe):
    \"\"\"
    Run Fisher exact test (overrepresentation) for every Biological Process
    GO term; return the minimum (most significant) BH-FDR across all terms.

    Returns np.nan when no gene in `genes` is annotated to any BP term.
    \"\"\"
    genes_set  = set(genes)
    N          = len(universe)
    n_selected = len(genes_set)
    if n_selected == 0:
        return np.nan

    pvals = []
    for annotated in go2gene.values():
        n_sig = len(set(annotated) & genes_set)
        if n_sig == 0:
            continue                          # term not hit — skip
        n_ann  = len(annotated)
        table  = [[n_sig,              n_selected - n_sig],
                  [n_ann - n_sig,  N - n_selected - n_ann + n_sig]]
        _, p   = fisher_exact(table, alternative="greater")
        pvals.append(p)

    if not pvals:
        return np.nan

    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    return float(np.min(fdr))

print("min_fdr_bp defined.")
print("Quick test (9-gene window vs PG10 universe):")
test_genes = mapping["IDpg10"].iloc[:9].tolist()
print(f"  min FDR = {min_fdr_bp(go2gene_pg10, test_genes, universe_pg10):.4f}")
"""
))

# ── SECTION 1 ────────────────────────────────────────────────────────────────
cells.append(md(
"""---
## 4. Background Comparison
### Does it matter which universe you use?

**The experiment:**
Consecutive BSU168 gene IDs (e.g. BSU00010–BSU00200) tend to cluster by
operon and metabolic pathway because the BSU168 genome is ordered by
chromosomal position.  A sliding window of 10 consecutive IDs is therefore
likely to hit a real functional group — making it a realistic proxy for an
actual gene list from a transcriptomics experiment.

**The two competing backgrounds:**

| Background | Universe | Justification |
|---|---|---|
| **PG10 correct** | 2,683 PG10 IDpg10 IDs | Matches the actual strain being studied |
| **BSU168 wrong** | 4,185 BSU168 IDbsub IDs | Would be used if the analyst did not account for the deletion strain |

**Expected outcome:**
Using a larger background (BSU168) makes each GO term appear less enriched
relative to the universe, which can either inflate FDR (false negatives)
or, because the universe denominator is different, also create false positives
when the wrong gene IDs are used.

We run **1,000 random sliding windows of 10 consecutive PG10-mapped BSU IDs**
and record the minimum FDR for each window under both backgrounds.
"""
))

cells.append(code(
"""# Sort the BSU IDs present in the PG10 mapping (2,678 unique IDs)
bsu_sorted  = sorted(mapping["IDbsub"].unique())
window_size = 10
n_windows   = len(bsu_sorted) - window_size + 1   # 2,669 possible windows

rng = np.random.default_rng(42)
starts = rng.choice(n_windows, size=1000, replace=False)

print(f"Unique BSU IDs in PG10 mapping : {len(bsu_sorted):,}")
print(f"Possible windows (size {window_size})    : {n_windows:,}")
print(f"Sampled windows                : {len(starts):,}")
"""
))

cells.append(code(
"""print("Running 1,000 windows × 2 backgrounds — this takes ~1–2 minutes...\\n")

fdr_pg10 = np.full(len(starts), np.nan)
fdr_bsub = np.full(len(starts), np.nan)

for i, start in enumerate(starts):
    bsu_win = bsu_sorted[start : start + window_size]

    # Correct: convert BSU IDs to PG10 IDpg10, test against PG10 universe
    pg10_genes   = mapping.loc[mapping["IDbsub"].isin(bsu_win), "IDpg10"].tolist()
    fdr_pg10[i]  = min_fdr_bp(go2gene_pg10, pg10_genes, universe_pg10)

    # Wrong: use BSU IDs directly against the full BSU168 universe
    bsu_genes    = [g for g in bsu_win if g in universe_bsub]
    fdr_bsub[i]  = min_fdr_bp(go2gene_bsub, bsu_genes, universe_bsub)

    if (i + 1) % 200 == 0:
        print(f"  {i+1:4d} / {len(starts)}")

# NA (no annotations in that window) → 1.0 (no enrichment)
fdr_pg10 = np.where(np.isnan(fdr_pg10), 1.0, fdr_pg10)
fdr_bsub = np.where(np.isnan(fdr_bsub), 1.0, fdr_bsub)

print("\\nDone.")
"""
))

cells.append(code(
"""# ── Quadrant summary ─────────────────────────────────────────────────────────
both     = ((fdr_pg10 < 0.05) & (fdr_bsub < 0.05)).sum()
pg10_only = ((fdr_pg10 < 0.05) & (fdr_bsub >= 0.05)).sum()
bsub_only = ((fdr_pg10 >= 0.05) & (fdr_bsub < 0.05)).sum()
neither  = ((fdr_pg10 >= 0.05) & (fdr_bsub >= 0.05)).sum()
n        = len(starts)

summary = pd.DataFrame({
    "Category"  : ["Both significant (FDR < 0.05)",
                   "PG10 only  → false negatives with BSU168",
                   "BSU168 only → false positives with BSU168",
                   "Neither significant"],
    "N windows" : [both, pg10_only, bsub_only, neither],
    "%" :         [round(v/n*100, 1) for v in [both, pg10_only, bsub_only, neither]],
})
print(summary.to_string(index=False))
print(f"\\nWindows with min-FDR < 0.20 :  PG10 = {(fdr_pg10 < 0.20).sum()}  |  BSU168 = {(fdr_bsub < 0.20).sum()}")
print(f"Windows with min-FDR = 1.00 :  PG10 = {(fdr_pg10 == 1).sum()}  |  BSU168 = {(fdr_bsub == 1).sum()}")
"""
))

cells.append(md(
"""### 4.1 R-generated figures (reference)

The figures below were produced by the original `background_comparison.R` script.
The Python reproductions follow.
"""
))

cells.append(code(
"""display(Image("background_effect_density.png", width=680))
"""
))

cells.append(code(
"""display(Image("background_effect_scatter.png", width=520))
"""
))

cells.append(md("### 4.2 Python reproductions"))

cells.append(code(
"""# ── Figure 1: overlapping density — FDR distributions ────────────────────────
col_pg10 = "#2166ac"   # blue  = correct
col_bsub = "#d73027"   # red   = wrong

fig, ax = plt.subplots(figsize=(8, 4.5))

for data, label, color in [
        (fdr_pg10, "PG10 correct  (n = 2,683)", col_pg10),
        (fdr_bsub, "BSU168 wrong  (n = 4,185)", col_bsub)]:
    ax.hist(data, bins=50, density=True, alpha=0.35,
            color=color, label=label, edgecolor="none")

    # KDE overlay
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(data, bw_method=0.15)
    xs  = np.linspace(0, 1, 300)
    ax.plot(xs, kde(xs), color=color, lw=2)

ax.axvline(0.05, ls="--", color="0.4", lw=1.2)
ax.text(0.067, ax.get_ylim()[1] * 0.95, "FDR 0.05",
        color="0.4", fontsize=9, va="top")
ax.set_xlim(0, 1)
ax.set_xlabel("Minimum FDR per window (lowest = most significant BP term)")
ax.set_ylabel("Density")
ax.set_title("Background set effect on GO Biological Process enrichment\\n"
             f"1,000 sliding windows of {window_size} consecutive PG10-mapped BSU168 genes")
ax.legend(loc="upper center", fontsize=10)
plt.tight_layout()
plt.savefig("background_effect_density_py.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# ── Figure 2: scatter — per-window PG10 FDR vs BSU168 FDR ───────────────────
quad_cols = {"Both":        "#4dac26",
             "PG10 only":   col_pg10,
             "BSU168 only": col_bsub,
             "Neither":     "0.7"}

quad = np.where((fdr_pg10 < 0.05) & (fdr_bsub < 0.05),  "Both",
       np.where((fdr_pg10 < 0.05) & (fdr_bsub >= 0.05), "PG10 only",
       np.where((fdr_pg10 >= 0.05) & (fdr_bsub < 0.05), "BSU168 only",
                                                          "Neither")))

fig, ax = plt.subplots(figsize=(5.5, 5.5))
for label, color in quad_cols.items():
    mask = quad == label
    ax.scatter(fdr_pg10[mask], fdr_bsub[mask],
               color=color, alpha=0.5, s=15, label=label, zorder=2)

ax.plot([0, 1], [0, 1], ls="--", color="0.3", lw=1, zorder=1)
ax.axhline(0.05, ls=":", color=col_bsub,  lw=1)
ax.axvline(0.05, ls=":", color=col_pg10,  lw=1)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel("Min FDR — PG10 background (correct)")
ax.set_ylabel("Min FDR — BSU168 background (wrong)")
ax.set_title("Per-window FDR: correct vs wrong background\\n"
             "Dashed diagonal = perfect agreement  |  dotted lines = FDR 0.05")
ax.legend(title="Significant at FDR 0.05", fontsize=9, loc="lower right")
plt.tight_layout()
plt.savefig("background_effect_scatter_py.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(md(
"""### 4.3 Interpretation

**Key takeaway:** Using the wrong background (BSU168 instead of PG10) causes
both false negatives and false positives.

- **False negatives (PG10-only significant):** The larger BSU168 universe
  contains genes absent from PG10, which dilutes the apparent enrichment
  signal.  A GO term that looks enriched in the PG10 context may not meet
  the FDR threshold when tested against a bigger background — even though
  the biological signal is real.

- **False positives (BSU168-only significant):** Using raw BSU IDs against
  the BSU168 universe bypasses the ID conversion step.  When a BSU168 gene
  has no PG10 counterpart, it should not appear in the universe at all.
  Its presence can accidentally enrich certain terms.

**Practical rule:** Always define the background universe as the set of genes
*that could have been detected in your experiment*, not the entire annotated
proteome of a reference strain.
"""
))

# ── SECTION 2 ────────────────────────────────────────────────────────────────
cells.append(md(
"""---
## 5. Sensitivity Analysis
### How much does a single gene change matter?

Small gene lists (< 20 genes) are common in targeted follow-up experiments
— knock-out sets, co-immunoprecipitation hits, operons of interest.  When
working near the FDR 0.05 boundary, removing or swapping just one gene can
completely rearrange the top-10 GO terms.

**Four scenarios** (all using the BSU168 universe for consistency with the R script):

| Set | Genes | Description |
|---|---|---|
| **A** | BSU00420–00500 (9 unique) | Full set |
| **B** | A minus BSU00430 (8 genes) | Drop one gene |
| **C** | A minus BSU00500 (8 genes) | Drop a different gene |
| **D** | A minus BSU00430 + BSU00510 (9 genes) | Swap one gene |

These BSU IDs correspond to consecutive genes in the BSU168 chromosome, all
annotated to overlapping GO Biological Process terms related to DNA
replication and repair.  Sets B and C each remove one gene; set D replaces
one with the next gene downstream in the same chromosomal region.

**Figure structure:** Two side-by-side heatmap panels share the same row axis
(the union of all GO terms appearing in any set's top-10 result):
- Left panel: **presence / absence** — is this term in the top-10 for this set?
- Right panel: **FDR value** (−log₁₀ scale) — how significant is it?
"""
))

cells.append(code(
"""def top_bp(genes, go2gene, universe, top_n=10):
    \"\"\"
    Return the top `top_n` GO Biological Process terms (by FDR ascending)
    for a given gene set.

    Returns None if no annotated gene is in the set.
    Returns a pd.DataFrame with columns: GO.ID, Term, pvalue, FDR.
    \"\"\"
    genes_set = set(g for g in genes if g in set(universe))
    N         = len(universe)
    n_s       = len(genes_set)
    if n_s == 0:
        return None

    rows = []
    for go_id, annotated in go2gene.items():
        n_sig = len(set(annotated) & genes_set)
        if n_sig == 0:
            continue
        n_ann = len(annotated)
        table = [[n_sig,           n_s - n_sig],
                 [n_ann - n_sig,  N - n_s - n_ann + n_sig]]
        _, p  = fisher_exact(table, alternative="greater")
        rows.append({"GO.ID": go_id, "pvalue": p})

    if not rows:
        return None

    res        = pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)
    _, fdr, _, _ = multipletests(res["pvalue"], method="fdr_bh")
    res["FDR"] = fdr
    res["Term"] = res["GO.ID"].map(go_names).fillna("")
    return res.head(top_n)

print("top_bp defined.")
"""
))

cells.append(code(
"""# ── Define the four gene sets ─────────────────────────────────────────────────
sets = {
    "A\\nfull\\n(n=9)"          : ["BSU00420","BSU00430","BSU00440","BSU00450",
                                    "BSU00460","BSU00470","BSU00480","BSU00490","BSU00500"],
    "B\\n−BSU00430\\n(n=8)"     : ["BSU00420","BSU00440","BSU00450","BSU00460",
                                    "BSU00470","BSU00480","BSU00490","BSU00500"],
    "C\\n−BSU00500\\n(n=8)"     : ["BSU00420","BSU00430","BSU00440","BSU00450",
                                    "BSU00460","BSU00470","BSU00480","BSU00490"],
    "D\\n−BSU00430\\n+BSU00510\\n(n=9)": ["BSU00420","BSU00440","BSU00450","BSU00460",
                                           "BSU00470","BSU00480","BSU00490","BSU00500","BSU00510"],
}

print("Running GO BP enrichment for 4 gene sets (BSU168 background)...\\n")
results = {nm: top_bp(genes, go2gene_bsub, universe_bsub, top_n=10)
           for nm, genes in sets.items()}

for nm, res in results.items():
    n = 0 if res is None else len(res)
    print(f"  {nm.replace(chr(10),' '):40s}: {n} terms returned")
"""
))

cells.append(code(
"""# ── Build union of all GO terms and long-format DataFrame ────────────────────
all_ids = []
seen    = set()
for res in results.values():
    if res is not None:
        for go_id in res["GO.ID"]:
            if go_id not in seen:
                all_ids.append(go_id)
                seen.add(go_id)

# Term label: "Name (GO:XXXXXXX)", truncated to 50 chars
def term_label(go_id):
    name = go_names.get(go_id, go_id)
    lbl  = f"{name} ({go_id})"
    return lbl[:47] + "..." if len(lbl) > 50 else lbl

# Row order: most frequent across sets first, then by median FDR
term_freq = {go_id: sum(r is not None and go_id in r["GO.ID"].values
                        for r in results.values())
             for go_id in all_ids}
term_mfdr = {}
for go_id in all_ids:
    fdrs = [r.loc[r["GO.ID"]==go_id,"FDR"].values[0]
            for r in results.values()
            if r is not None and go_id in r["GO.ID"].values]
    term_mfdr[go_id] = np.median(fdrs) if fdrs else 1.0

ordered_ids = sorted(all_ids,
    key=lambda x: (-term_freq[x], term_mfdr[x]))
ordered_labels = [term_label(g) for g in ordered_ids]

# Long-format DataFrame
rows = []
set_names = list(sets.keys())
for nm in set_names:
    res = results[nm]
    for go_id, lbl in zip(ordered_ids, ordered_labels):
        if res is not None and go_id in res["GO.ID"].values:
            fdr = res.loc[res["GO.ID"]==go_id,"FDR"].values[0]
            rows.append({"term": lbl, "set": nm, "present": True,  "fdr": fdr})
        else:
            rows.append({"term": lbl, "set": nm, "present": False, "fdr": np.nan})

df_long = pd.DataFrame(rows)
df_long["log10fdr"] = -np.log10(df_long["fdr"])
df_long["term"]     = pd.Categorical(df_long["term"], categories=list(reversed(ordered_labels)))
df_long["set"]      = pd.Categorical(df_long["set"],  categories=set_names)

print(f"Union of GO terms: {len(ordered_ids)}")
print(f"Long-format rows : {len(df_long)}")
"""
))

cells.append(md(
"""### 5.1 R-generated figure (reference)

The figure below was produced by `sensitivity_analysis.R`.
The Python reproduction follows.
"""
))

cells.append(code(
"""display(Image("sensitivity_heatmap.png", width=780))
"""
))

cells.append(md("### 5.2 Python reproduction"))

cells.append(code(
"""n_rows  = len(ordered_ids)
fig_h   = max(4.5, 1.3 + n_rows * 0.35)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, fig_h),
                                gridspec_kw={"width_ratios": [1.6, 1]})

# ── Panel 1: presence / absence ──────────────────────────────────────────────
piv1 = df_long.pivot(index="term", columns="set", values="present")
piv1 = piv1.loc[list(reversed(ordered_labels))]   # rows: top terms first

im1 = ax1.imshow(
    piv1.values.astype(float),
    aspect="auto", cmap="Blues", vmin=0, vmax=1,
    extent=[-0.5, len(set_names)-0.5, -0.5, n_rows-0.5]
)
ax1.set_xticks(range(len(set_names)))
ax1.set_xticklabels(
    [nm.replace("\\n", " ") for nm in set_names],
    fontsize=8, ha="center"
)
ax1.set_yticks(range(n_rows))
ax1.set_yticklabels(list(reversed(ordered_labels)), fontsize=7.5)
ax1.xaxis.set_label_position("top"); ax1.xaxis.tick_top()

# grid lines
for x in np.arange(-0.5, len(set_names), 1):
    ax1.axvline(x, color="white", lw=1)
for y in np.arange(-0.5, n_rows, 1):
    ax1.axhline(y, color="white", lw=0.5)

ax1.set_title("Presence / absence\\n(in top 10)", fontsize=10, fontweight="bold", pad=8)

# legend patches
ax1.legend(
    handles=[mpatches.Patch(color="#2166ac", label="in top 10"),
             mpatches.Patch(color="#deebf7", label="absent")],
    loc="lower left", fontsize=8, framealpha=0.8
)

# ── Panel 2: FDR value ────────────────────────────────────────────────────────
piv2 = df_long.pivot(index="term", columns="set", values="log10fdr")
piv2 = piv2.loc[list(reversed(ordered_labels))]

# Build masked array so NaN shows as grey
import numpy.ma as ma
masked = ma.masked_invalid(piv2.values.astype(float))

cmap = plt.cm.Blues.copy()
cmap.set_bad("lightgrey")

im2 = ax2.imshow(masked, aspect="auto", cmap=cmap,
                 vmin=0, vmax=np.nanmax(piv2.values) or 1,
                 extent=[-0.5, len(set_names)-0.5, -0.5, n_rows-0.5])
ax2.set_xticks(range(len(set_names)))
ax2.set_xticklabels(["A","B","C","D"], fontsize=10, fontweight="bold")
ax2.set_yticks([])
ax2.xaxis.set_label_position("top"); ax2.xaxis.tick_top()

for x in np.arange(-0.5, len(set_names), 1):
    ax2.axvline(x, color="white", lw=1)
for y in np.arange(-0.5, n_rows, 1):
    ax2.axhline(y, color="white", lw=0.5)

ax2.set_title("FDR value\\n(dark = significant)", fontsize=10, fontweight="bold", pad=8)

cbar = fig.colorbar(im2, ax=ax2, orientation="horizontal",
                    pad=0.04, fraction=0.04, aspect=20)
cbar.set_label(r"$-\log_{10}$(FDR)", fontsize=9)

# ── Suptitle ─────────────────────────────────────────────────────────────────
fig.suptitle(
    "Sensitivity of GO BP enrichment to single-gene changes\\n"
    f"BSU168 background (n = {len(universe_bsub):,})  ·  "
    f"top 10 terms per set  ·  {len(all_ids)} unique terms shown\\n"
    "Set A: " + ", ".join(sets[list(sets.keys())[0]]),
    fontsize=10, y=1.01
)
plt.tight_layout()
plt.savefig("sensitivity_heatmap_py.png", dpi=150, bbox_inches="tight")
plt.show()
"""
))

cells.append(code(
"""# ── Save results table ────────────────────────────────────────────────────────
out_rows = []
for nm, res in results.items():
    if res is not None:
        r = res.copy()
        r["set"] = nm.replace("\\n", " ")
        out_rows.append(r[["set","GO.ID","Term","pvalue","FDR"]])

df_out = pd.concat(out_rows, ignore_index=True)
df_out.to_csv("sensitivity_analysis_results_py.csv", index=False)
print(f"Saved: sensitivity_analysis_results_py.csv  ({len(df_out)} rows)")
print()
print(df_out.to_string(index=False))
"""
))

cells.append(md(
"""### 5.3 Interpretation

**Key takeaway:** Near the FDR 0.05 boundary, GO enrichment results are
highly sensitive to the exact gene set composition.

- **Set B (−BSU00430):** Removing a single gene drops several terms that were
  significant in Set A, while others become newly significant or shift rank.
- **Set C (−BSU00500):** A different gene removal produces a different
  rearrangement — the affected terms differ from Set B.
- **Set D (swap BSU00430 → BSU00510):** The replacement gene lies in the
  same chromosomal region and carries overlapping GO annotations, yet the
  resulting term landscape still differs from Set A.

This sensitivity arises because Fisher exact tests with small gene sets
(here n ≤ 9) have very limited statistical power.  The BH correction
distributes any p-value shift across all tested terms simultaneously, so
a single gene that "bridges" two GO terms can flip several FDRs across the
0.05 threshold at once.

**Practical implications:**
1. Report enrichment results with confidence intervals or permutation-based
   p-values when the gene set is small (< 20 genes).
2. Treat terms that appear only marginally significant (0.01 < FDR < 0.10)
   as hypotheses requiring independent validation, not conclusions.
3. Consider testing multiple background sets and reporting only terms
   consistent across reasonable background choices.
"""
))

cells.append(md(
"""---
## 6. Summary

| Analysis | Core finding |
|---|---|
| **Background comparison** | Using the wrong (larger) universe creates false negatives *and* false positives; the PG10-matched background is essential for PG10 experiments |
| **Sensitivity analysis** | Single-gene changes near FDR 0.05 can completely reshuffle the top-10 GO terms; small-set enrichment results should be treated as exploratory |

**Technical note:** Both analyses implement GO enrichment from scratch using
`scipy.stats.fisher_exact` + `statsmodels.stats.multitest.multipletests`
(Benjamini–Hochberg).  This is equivalent to topGO's classic/Fisher algorithm
and to the `runGOFisher` function in `app.R`, making the statistical machinery
fully transparent and reproducible without any GO-specific library.
"""
))

# ── BUILD NOTEBOOK ───────────────────────────────────────────────────────────
nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {
            "display_name": "Python (plant_immunity)",
            "language": "python",
            "name": "plant_immunity"
        },
        "language_info": {
            "name": "python",
            "version": "3.12.0"
        }
    },
    "cells": cells,
}

out = pathlib.Path(__file__).with_name("BsubtGO_analysis.ipynb")
with open(out, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f"Notebook written: {out}")
print(f"Total cells: {len(cells)}  "
      f"(md: {sum(1 for c in cells if c['cell_type']=='markdown')}  "
      f"code: {sum(1 for c in cells if c['cell_type']=='code')})")
