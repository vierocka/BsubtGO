# BsubtilisGO

A Shiny web application for Gene Ontology (GO) enrichment analysis of *Bacillus subtilis* gene sets, designed for wet-lab biologists who work with strains other than the standard reference.

Live app: **https://vierakovacova.shinyapps.io/BsubtilisGO/**

---

## The problem this solves

GO enrichment analysis asks: are any biological functions statistically over-represented in my gene list compared to what I would expect by chance? The answer depends critically on what "by chance" means — i.e., on the **background gene set** used as the reference universe.

Most online GO tools (g:Profiler, DAVID, STRING web interface) only carry *B. subtilis* strain **168** in their databases. If your experiment was done in a different strain, using the 168 background introduces two errors:

1. **Genes present in 168 but absent from your strain** are included in the background, diluting enrichment signals and producing false negatives.
2. **Genes specific to your strain with no 168 equivalent** cannot be entered at all, or are silently dropped, distorting the statistics.

This is particularly relevant for engineered or reduced-genome strains whose protein-coding repertoire differs meaningfully from strain 168.

**BsubtilisGO** addresses this by building a strain-matched background: GO annotations are transferred from the well-annotated strain 168 to each target strain via BLASTP best-hit homology, and the enrichment test is run against only the genes that actually exist in the strain of interest.

---

## Supported strains

| Strain | ID format | Background size |
|--------|-----------|-----------------|
| BSUB168 | `BSU*` (e.g. `BSU00240`) | 4,185 proteins |
| PG10 | `ANY*.1` (e.g. `ANY33920.1`) | 2,769 proteins |

---

## Usage

The app is intentionally minimal so it is quick to use at the bench.

1. Select your strain from the dropdown.
2. Paste your gene IDs — comma- or newline-separated — into the text box.  
   The placeholder updates automatically to show the correct ID format for the selected strain.
3. Click **Analyze**.  
   Results appear as three sortable tables: Biological Process, Molecular Function, and Cellular Component.
4. Click **Download Results** to export all three tables as a single CSV.

**Example IDs:**

| Strain  | Example input |
|---------|---------------|
| BSUB168 | `BSU00240, BSU00260, BSU00280, BSU00290, BSU00300` |
| PG10    | `ANY33920.1, ANY33921.1, ANY33922.1, ANY33923.1` |

---

## Methods

### 1. Reference proteome and GO annotations — *B. subtilis* 168

The proteome of *B. subtilis* 168 (STRING-db species ID 224308) was downloaded from STRING-db v12.0:

| File | Contents |
|------|----------|
| `224308.protein.sequences.v12.0.fa` | 4,185 protein sequences (FASTA) |
| `224308.protein.enrichment.terms.v12.0.txt` | Functional annotations per protein: GO terms, UniProt keywords, Pfam/InterPro domains, Reactome pathways, subcellular localisation (COMPARTMENTS) |
| `224308.protein.aliases.v12.0.txt` | Cross-reference aliases (RefSeq, UniProt, etc.) |
| `224308.clusters.proteins.v12.0.txt` | STRING cluster assignments |

STRING protein IDs use the format `224308.BSU00010` (taxonomy ID + BSU locus tag). Only the three GO categories were used in this project:

- `Biological Process (Gene Ontology)`
- `Molecular Function (Gene Ontology)`
- `Cellular Component (Gene Ontology)`

---

### 2. *B. subtilis* PG10 proteome

PG10 protein sequences and genome annotations were retrieved from NCBI (accession **CP016788**):

| File | Contents |
|------|----------|
| `NCBI_PG10_sequences_aa.fa` | 2,770 protein sequences (FASTA); headers contain NCBI protein IDs (e.g. `ANY33920.1`) |
| `Bacillus_subtilis_PG10_reducedGenome.gff3` | Genome annotation in GFF3 format |
| `Proteins_inPG10.gtf` | Protein features in GTF format |
| `Proteins_inPG10_referenceIDs.list` | WP_* reference protein IDs for PG10 genes |

Protein names were extracted from FASTA headers using `protein_names_list_make.sh`:

```bash
grep ">" NCBI_PG10_sequences_aa.fa \
  | awk 'BEGIN{FS="["}; { print $3}' \
  | sed s/protein=// \
  | tr -d '\]' > protein_names.list
```

---

### 3. BLASTP — mapping PG10 proteins to BSU168 homologs

A protein BLAST database was built from the BSU168 reference proteome:

```bash
makeblastdb -in 224308.protein.sequences.v12.0.fa -dbtype prot -out 224308.protein.db
```

PG10 proteins were aligned against this database to identify the single best-hit BSU168 homolog for each PG10 gene:

```bash
blastp \
  -db    224308.protein.db \
  -query NCBI_PG10_sequences_aa.fa \
  -outfmt 6 \
  -max_target_seqs 1 \
  -max_hsps 1 \
  -out   PG10ncbiProt_stringBsubtRefDB.blp
```

Output format: BLAST tabular format 6 (`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`).

**Results summary:** 2,769 hits for 2,770 query proteins. The majority are near-identical (most e-values = 0, percent identity = 100 %), consistent with the close relationship between PG10 and strain 168. A small number of hits have low sequence identity (minimum ~18.6 %) and poor e-values (up to ~8), likely representing PG10-specific genes with no true BSU168 ortholog. **No e-value cutoff was applied** — all best hits were retained. GO terms inherited by poor-match entries should be interpreted with caution.

---

### 4. Building the GO mapping table

The shell script `pipeline_ID_conversion.txt` joins the BLAST output with the STRING-db enrichment terms to produce the mapping table used by the app. For each BLAST hit it:

1. Extracts the PG10 protein ID from the query field  
   (e.g. `lcl|CP016788.1_prot_ANY33920.1_1` → `ANY33920.1`)
2. Extracts the BSU168 locus tag from the subject field  
   (e.g. `224308.BSU00010` → `BSU00010`)
3. Looks up all GO terms in `224308.protein.enrichment.terms.v12.0.txt` for each of the three categories using `awk`
4. Concatenates multiple GO terms with a space separator

```bash
echo -e "IDpg10\tIDbsub\tBiolProc\tMolFun\tCellComp" \
  > PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv

count=$(cut -f1 PG10ncbiProt_stringBsubtRefDB.blp | sort | wc -l | cut -d" " -f1)

for ((i=1; i<=$(($count)); i++)); do
  IDpg10=$(sed -n "${i}p" PG10ncbiProt_stringBsubtRefDB.blp | cut -f1 | cut -d"_" -f3)
  IDbsub=$(sed -n "${i}p" PG10ncbiProt_stringBsubtRefDB.blp | cut -f2 | cut -d"." -f2)
  BiolP=$(awk  'BEGIN{FS="\t"}; /'"$IDbsub"'/ && /Biological Process/ { print $3}' \
           224308.protein.enrichment.terms.v12.0.txt)
  MolF=$(awk   'BEGIN{FS="\t"}; /'"$IDbsub"'/ && /Molecular Function/ { print $3}' \
           224308.protein.enrichment.terms.v12.0.txt)
  Compart=$(awk 'BEGIN{FS="\t"}; /'"$IDbsub"'/ && /Cellular Component/ { print $3}' \
           224308.protein.enrichment.terms.v12.0.txt)
  echo -e "$IDpg10\t$IDbsub\t$BiolP\t$MolF\t$Compart" \
    >> PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv
done
```

Output: `PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv` — 2,769 gene rows, five tab-separated columns:

| Column | Description |
|--------|-------------|
| `IDpg10` | PG10 NCBI protein ID (e.g. `ANY33920.1`) |
| `IDbsub` | BSU168 locus tag (e.g. `BSU00010`) |
| `BiolProc` | Space-separated GO terms — Biological Process |
| `MolFun` | Space-separated GO terms — Molecular Function |
| `CellComp` | Space-separated GO terms — Cellular Component |

Genes with no annotation in a given category have an empty field (877 in BiolProc, 802 in MolFun, 1,141 in CellComp). The app handles these silently.

> **Note:** Regenerating this file with DuckDB is much faster than the bash loop above — DuckDB can join `.blp` and `enrichment.terms` in a single SQL query in under a second, versus several minutes for the sequential `awk` approach.

---

### 5. GO enrichment analysis — Shiny app (`app.R`)

The app implements GO over-representation analysis using base R only (`shiny` and `DT` are the sole dependencies, both from CRAN). No Bioconductor packages are required.

On startup the app reads the mapping CSV and pre-computes gene-to-GO list objects for both strains. GO term descriptions are loaded from `go_terms.tsv`, a 4,065-term lookup table extracted from `224308.protein.enrichment.terms.v12.0.txt`.

For each GO term, a one-sided Fisher exact test is run via `fisher.test(alternative = "greater")`, testing whether the term is over-represented in the submitted gene list relative to the strain background. This is equivalent to the classic/Fisher algorithm in `topGO`.

| Parameter | Value |
|-----------|-------|
| Test | Fisher exact test, one-sided (over-representation) |
| Background | All annotated genes in the selected strain |
| Multiple testing correction | FDR (Benjamini–Hochberg), applied to the top 20 results per ontology |
| Terms reported | Top 20 per ontology (by p-value) |
| Strains | BSUB168 (`BSU*`) and PG10 (`ANY*.1`) |
| Dependencies | `shiny`, `DT` (CRAN only) |

> **Limitation:** FDR is calculated on the truncated top-20 list rather than across all tested GO terms. It is indicative but not a genome-wide multiple-testing correction. Treat borderline FDR values with caution.

---

### 6. Auxiliary — KEGG pathway analysis (`Kegg.R`)

A standalone exploratory script (not part of the Shiny app) that queries the KEGG database for *B. subtilis* pathway enrichment using `clusterProfiler::enrichKEGG` (organism code `bsu`) and visualises selected pathways with `pathview`.

---

## Deployment

The app depends only on CRAN packages (`shiny`, `DT`) — no Bioconductor required.

```r
library(rsconnect)
rsconnect::setAccountInfo(
  name   = "vierakovacova",
  token  = "YOUR_TOKEN",   # shinyapps.io → Account → Tokens
  secret = "YOUR_SECRET"
)
rsconnect::deployApp(
  appDir   = "path/to/BsubtGO",
  appName  = "BsubtilisGO",
  appFiles = c("app.R",
               "PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv",
               "go_terms.tsv",
               ".Rprofile")
)
```

---

## File inventory

```
BsubtGO/                                          ← this repository
├── app.R                                         ← Shiny app
├── PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv  ← GO mapping table (deploy with app)
├── go_terms.tsv                                  ← GO ID → term name lookup (deploy with app)
├── .Rprofile                                     ← CRAN repo config
└── README.md

help_files/                                       ← source data and pipeline (not deployed)
├── 224308.protein.sequences.v12.0.fa             ← BSU168 proteome (STRING-db v12.0)
├── 224308.protein.enrichment.terms.v12.0.txt     ← BSU168 GO + functional annotations
├── 224308.protein.aliases.v12.0.txt              ← BSU168 cross-reference aliases
├── 224308.clusters.proteins.v12.0.txt            ← STRING cluster assignments
├── 224308.protein.db.*                           ← BLAST protein DB (built from BSU168 FASTA)
├── NCBI_PG10_sequences_aa.fa                     ← PG10 proteome (NCBI CP016788)
├── Bacillus_subtilis_PG10_reducedGenome.gff3     ← PG10 genome annotation
├── Proteins_inPG10.gtf                           ← PG10 protein features (GTF)
├── Proteins_inPG10_referenceIDs.list             ← PG10 WP_* reference IDs
├── protein_names.list                            ← PG10 protein names (from FASTA headers)
├── protein_list.txt                              ← PG10 protein ID list
├── protein_names_list_make.sh                    ← extracts protein names from FASTA headers
├── PG10ncbiProt_stringBsubtRefDB.blp             ← BLASTP output (format 6)
├── pipeline_ID_conversion.txt                    ← bash pipeline that builds the mapping CSV
├── PG10_flgM_mScarlet_Oct2024.bed.mat            ← genomics experiment data (Oct 2024)
└── Kegg.R                                        ← standalone KEGG pathway analysis script
```
