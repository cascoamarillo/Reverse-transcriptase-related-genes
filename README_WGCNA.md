# WGCNA co-expression analysis - *Adineta ricciae*

Weighted gene co-expression network analysis (WGCNA) of the bdelloid rotifer
*Adineta ricciae*, re-analysing the RNA-seq data of **Nowell et al. 2024**
(*Nat. Commun.* 15:5787, https://doi.org/10.1038/s41467-024-49919-1). The raw
reads and the standard DESeq2 differential-expression analysis are from that
study; this repository adds the WGCNA module-level analysis.

## Experimental design

12 libraries: fungal-infection time course, 4 groups × 3 biological replicates.

| Group | Infection | Time (hpi) | Samples |
|-------|-----------|-----------|---------|
| Control_7h  | control  | 7  | ctrl.7a, ctrl.7b, ctrl.7c |
| Control_24h | control  | 24 | ctrl.24b, ctrl.24c, ctrl.24d |
| Infected_7h  | infected | 7  | treatment.7a, treatment.7c, treatment.7d |
| Infected_24h | infected | 24 | treatment.24b, treatment.24c, treatment.24d |

Sample name pattern: `<ctrl|treatment>.<7|24><replicate letter>`. Gene IDs are `ARIC|g#`.

## Workflow

The pipeline (`wgcna_Ar.R`) follows the ISU Bioinformatics Workbook WGCNA
tutorial (https://bioinformaticsworkbook.org/tutorials/wgcna.html), adapted to
the *A. ricciae* count matrix:

1. **Load** the gene-level count matrix (genes × samples).
2. **Normalize** with DESeq2 (`DESeqDataSetFromMatrix` → `DESeq`), design `~ Group`.
3. **Transform / filter** — variance-stabilizing transformation (`getVarianceStabilizedData`)
   and retain the top-variance genes (default: 95th-percentile variance, as in the tutorial).
4. **Soft threshold** — `pickSoftThreshold` over powers 1–40; pick the lowest
   power reaching signed scale-free R² ≥ 0.90.
5. **Network / modules** — `blockwiseModules` (signed network, `minModuleSize = 30`,
   `mergeCutHeight = 0.25`, `deepSplit = 2`).
6. **Module colors + dendrogram** (`labels2colors`, `plotDendroAndColors`) → `Ar_gene_modules.txt`.
7. **Module eigengenes** (`moduleEigengenes` → `orderMEs`) — a module-by-sample
   heatmap plus a **module-trait correlation** heatmap (`labeledHeatmap`) against
   infection status, timepoint, and late (24 hpi) infection → `Ar_module_trait_cor.txt`.
8. **Modules of interest** — expression profiles and a topological-overlap
   (`TOMsimilarityFromExpr`) **edge list** for Cytoscape/VisANT → `Ar_edgelist.tsv`.

An optional appendix reproduces the alternative **manual** route
(`adjacency` → `TOMdist` → `cutreeDynamic`) on the DE-significant gene subset.

## Key parameters (set in the CONFIG block)

| Parameter | Default | Notes |
|-----------|---------|-------|
| `USE_VST_FILTER` | `FALSE` | `FALSE` reproduces the **published** figure (raw DESeq-normalized counts, all genes). `TRUE` = tutorial default (VST + variance filter), the more conventional WGCNA input. |
| `VAR_QUANTILE` | `0.95` | keep genes above this variance quantile (used only if `USE_VST_FILTER = TRUE`) |
| `PICKED_POWER` | `20` | **published value.** Signed network, 12 samples (WGCNA recommends ~18 for <20 samples); scale-free R² ≈ 0.92, on the plateau. Turquoise module contains the *A. ricciae* RVT homolog. |
| `NETWORK_TYPE` | `signed` | |
| `MIN_MODULE_SIZE` | `30` | |
| `MERGE_CUTHEIGHT` | `0.25` | |
| `MODULES_OF_INTEREST` | `turquoise, blue` | turquoise = module containing the RVT homolog |

## Inputs

- `Ar_counts.1.txt` — gene-level counts (htseq-count), tab-delimited; first column
  `Count` = gene ID (`ARIC|g#`), then the 12 sample columns. The sample→group
  mapping is parsed automatically from the `<ctrl|treatment>.<7|24><letter>` names
  in section 2 of the script (no editing needed for this dataset).

## Outputs

- `Ar_normalized_counts.txt` — DESeq2-normalized counts
- `Ar_gene_modules.txt` — gene → module-color assignments
- `Ar_module_trait_cor.txt` — module-eigengene correlations with infection/time traits
- `Ar_TOM-block.*.RData` — saved topological-overlap matrices
- `Ar_edgelist.tsv` — edge list (modules of interest) for Cytoscape/VisANT
- soft-threshold, dendrogram, and module-eigengene plots

## Deviations from the tutorial (documented)

- The tutorial fed a **VST + 95%-variance-filtered** matrix at **power 9**
  (maize, 24 samples). Here a **signed** network on **12 samples** uses
  **power 20** (WGCNA recommends ~18 for <20 samples; scale-free R² ≈ 0.92, on
  the plateau), which resolves differentiated modules — the **turquoise** module
  contains the *A. ricciae* RVT homolog. The published figure used raw
  DESeq-normalized counts (`USE_VST_FILTER = FALSE`); flip to `TRUE` for the
  VST/variance-filtered variant.
- Duplicated runs, leftover tutorial example output, and interactive `?help()`
  calls from the original working script have been removed.

## Reproducing

```r
# edit the CONFIG block at the top of wgcna_Ar.R, then:
source("wgcna_Ar.R")
```

## Software

R with `tidyverse`, `magrittr`, `DESeq2`, `WGCNA`, `genefilter`
(record the exact package versions with `sessionInfo()` for the paper).

## Citations

- Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation
  network analysis. *BMC Bioinformatics* 9:559. https://doi.org/10.1186/1471-2105-9-559
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biol.* 15:550.
  https://doi.org/10.1186/s13059-014-0550-8
- Nowell RW, Rodriguez F et al. (2024). Bdelloid rotifers deploy horizontally acquired
  biosynthetic genes against a fungal pathogen. *Nat. Commun.* 15:5787.
  https://doi.org/10.1038/s41467-024-49919-1
- Tutorial: ISU Bioinformatics Workbook — WGCNA.
  https://bioinformaticsworkbook.org/tutorials/wgcna.html

## Data availability

Raw reads and DESeq2 DGE: see Nowell et al. 2024 (deposited under NCBI
BioProject RJEB39927 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB39927).
