# BsubtGO

A Shiny web application for Gene Ontology (GO) enrichment analysis of *Bacillus subtilis* gene sets. The app supports two strains — the reference strain **168** (NCBI taxonomy 224308) and the engineered strain **PG10** — and performs Fisher exact tests across all three GO ontologies (Biological Process, Molecular Function, Cellular Component) using the `topGO` R package.

---

## Usage

1. Select the strain whose gene IDs you will enter (BSUB168 or PG10).
2. Paste gene IDs separated by commas or newlines into the text box.
3. Click **Analyze**. Results for BP, MF, and CC appear as sortable tables.
4. Download the combined results as a CSV with **Download Results**.

**Example IDs:**

| Strain   | Example input                              |
|----------|--------------------------------------------|
| BSUB168  | `BSU00240, BSU00260, BSU00280`             |
| PG10     | `ANY33920.1, ANY33921.1, ANY33922.1`       |

---

## Methods

### 1. Reference proteome and GO annotations — *B. subtilis* 168

The proteome of *B. subtilis* 168 (STRING-db species ID 224308) was downloaded from STRING-db v12.0:

| File | Contents |
|------|----------|
| `224308.protein.sequences.v12.0.fa` | 4,185 protein sequences (FASTA) |
| `224308.protein.enrichment.terms.v12.0.txt` | Functional annotations per protein, including GO terms, UniProt keywords, Pfam/InterPro domains, Reactome pathways, and subcellular localisation (COMPARTMENTS) |
| `224308.protein.aliases.v12.0.txt` | Cross-reference aliases from RefSeq, UniProt, etc. |
| `224308.clusters.proteins.v12.0.txt` | STRING cluster assignments |

STRING protein IDs use the format `224308.BSU00010`, where the numeric prefix is the NCBI taxonomy ID and the suffix is the BSU locus tag. Only the three GO categories were used in this project:

- `Biological Process (Gene Ontology)`
- `Molecular Function (Gene Ontology)`
- `Cellular Component (Gene Ontology)`

---

### 2. *B. subtilis* PG10 proteome

PG10 protein sequences and genome annotations were retrieved from NCBI (accession **CP016788**):

| File | Contents |
|------|----------|
| `NCBI_PG10_sequences_aa.fa` | 2,770 protein sequences (FASTA); headers contain NCBI protein IDs (e.g., `ANY33920.1`) |
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

PG10 proteins were then aligned against this database to identify best-hit BSU168 homologs:

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

**Results summary (`PG10ncbiProt_stringBsubtRefDB.blp`):** 2,769 hits for 2,770 query proteins. The majority of hits are near-identical (most e-values = 0, percent identity = 100 %), reflecting the close relationship between PG10 and strain 168. A small number of hits have low sequence identity (minimum ~18.6 %) and poor e-values (up to ~8), likely corresponding to PG10-specific genes that lack a true BSU168 ortholog. **No e-value cutoff was applied**; all best hits were retained. GO terms inherited by these poor-match entries should be interpreted with caution.

---

### 4. Building the GO mapping table

The shell script `pipeline_ID_conversion.txt` joins the BLAST output with the STRING-db enrichment terms to produce the final mapping table used by the Shiny app. For each BLAST hit, it:

1. Extracts the PG10 protein ID from the query field  
   (e.g., `lcl|CP016788.1_prot_ANY33920.1_1` → `ANY33920.1`)
2. Extracts the BSU168 locus tag from the subject field  
   (e.g., `224308.BSU00010` → `BSU00010`)
3. Looks up all matching GO terms in `224308.protein.enrichment.terms.v12.0.txt` for each of the three GO categories, using `awk` pattern matching on the locus tag and category name
4. Concatenates multiple GO terms for the same gene with a space separator

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

Output: `PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv` — 2,769 gene rows, tab-separated, five columns:

| Column | Description |
|--------|-------------|
| `IDpg10` | PG10 NCBI protein ID (e.g., `ANY33920.1`) |
| `IDbsub` | BSU168 locus tag (e.g., `BSU00010`) |
| `BiolProc` | Space-separated GO terms — Biological Process |
| `MolFun` | Space-separated GO terms — Molecular Function |
| `CellComp` | Space-separated GO terms — Cellular Component |

Genes with no annotation in a given category have an empty field (877 in BiolProc, 802 in MolFun, 1,141 in CellComp). The Shiny app handles these gracefully.

> **Note:** Regenerating this file via DuckDB is much faster than the bash loop above, as DuckDB can join `.blp` and `enrichment.terms` directly in a single SQL query. The bash pipeline takes several minutes; the equivalent DuckDB query runs in under a second.

---

### 5. GO enrichment analysis — Shiny app (`app.R`)

The app uses the `topGO` R package. On startup it reads the mapping CSV and pre-computes gene-to-GO list objects for both strains so that per-user analysis is fast.

| Parameter | Value |
|-----------|-------|
| Test | Fisher exact test |
| Algorithm | classic (no topology correction) |
| Correction | FDR (Benjamini–Hochberg), applied to the top results returned by `GenTable` |
| Nodes reported | Top 20 per ontology |
| Strains | BSUB168 (`BSU*` IDs) and PG10 (`ANY*.1` IDs) |

The FDR is calculated on the truncated top-20 list rather than all tested GO terms, so it should be treated as indicative rather than a genome-wide multiple-testing correction.

---

### 6. Auxiliary — KEGG pathway analysis (`Kegg.R`)

A standalone exploratory script that queries the KEGG database for *B. subtilis* pathway enrichment using `clusterProfiler::enrichKEGG` (organism code `bsu`) and visualises selected pathways with `pathview`. This script is not part of the Shiny app.

---

## Deployment to shinyapps.io

`topGO` is a Bioconductor package. The `.Rprofile` in this directory sets Bioconductor repositories so shinyapps.io can install all dependencies automatically.

```r
# Run once inside the BsubtGO/ directory
install.packages(c("BiocManager", "rsconnect"))
rsconnect::deployApp()
```

---

## File inventory

```
BsubtGO/                                         ← deploy from here
├── app.R                                        ← Shiny app
├── PG10id_BSUBid_goBiolP_goMolF_goCellComp.csv ← GO mapping table
├── .Rprofile                                    ← Bioconductor repo config
├── README.md
└── rsconnect/                                   ← deployment metadata

help_files/                                      ← source data and pipeline
├── 224308.protein.sequences.v12.0.fa            ← BSU168 proteome (STRING-db v12.0)
├── 224308.protein.enrichment.terms.v12.0.txt    ← BSU168 GO + functional annotations
├── 224308.protein.aliases.v12.0.txt             ← BSU168 cross-reference aliases
├── 224308.clusters.proteins.v12.0.txt           ← STRING cluster assignments
├── 224308.protein.db.*                          ← BLAST protein DB (built from above)
├── NCBI_PG10_sequences_aa.fa                    ← PG10 proteome (NCBI CP016788)
├── Bacillus_subtilis_PG10_reducedGenome.gff3    ← PG10 genome annotation
├── Proteins_inPG10.gtf                          ← PG10 protein features (GTF)
├── Proteins_inPG10_referenceIDs.list            ← PG10 WP_* reference IDs
├── protein_names.list                           ← PG10 protein names extracted from FASTA
├── protein_list.txt                             ← PG10 protein ID list
├── protein_names_list_make.sh                   ← script to extract names from FASTA headers
├── PG10ncbiProt_stringBsubtRefDB.blp            ← BLASTP output (format 6)
├── pipeline_ID_conversion.txt                   ← bash pipeline → builds the mapping CSV
├── PG10_flgM_mScarlet_Oct2024.bed.mat           ← ChIP/genomics experiment data (Oct 2024)
└── Kegg.R                                       ← standalone KEGG pathway analysis script
```
