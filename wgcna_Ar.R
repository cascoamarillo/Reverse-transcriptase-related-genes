################################################################################
# WGCNA co-expression network analysis of the bdelloid rotifer Adineta ricciae
#
# Re-analysis of the RNA-seq data from Nowell, Rodriguez et al. 2024
#   (Nat. Commun. 15:5787; https://doi.org/10.1038/s41467-024-49919-1)
# adding weighted gene co-expression network analysis (WGCNA) to the
# standard DESeq2 differential-expression analysis.
#
# Workflow adapted from the ISU Bioinformatics Workbook WGCNA tutorial:
#   https://bioinformaticsworkbook.org/tutorials/wgcna.html
# Method: Langfelder & Horvath 2008 (https://doi.org/10.1186/1471-2105-9-559)
#
# Input : gene-level count matrix (genes x samples) for A. ricciae
# Output: soft-threshold plots, module assignments, module-eigengene heatmap,
#         and a TOM edge list for Cytoscape/VisANT.
################################################################################

## ============================================================================
## 0. CONFIG  — edit these, then run top to bottom
## ============================================================================
WORKDIR        <- "."                      # set to your working folder
COUNTS_FILE    <- "Ar_counts.txt"        # genes x samples, tab-delimited; col 1 = gene ID
OUT_PREFIX     <- "Ar"                     # prefix for output files

USE_VST_FILTER <- FALSE  # FALSE reproduces the PUBLISHED figure: raw DESeq-normalized counts,
                         #        all genes (the original power-20 run). TRUE = tutorial default
                         #        (VST + variance filter), the more conventional WGCNA input.
VAR_QUANTILE   <- 0.95   # keep genes above this variance quantile (used only if USE_VST_FILTER=TRUE)

PICKED_POWER   <- 20     # Published value. Signed network with 12 samples: WGCNA guidance
                         #        recommends power ~18 for <20 samples, and the scale-free fit is on
                         #        its R^2 >= 0.90 plateau here (R^2 ~ 0.92 at power 20). Power 20
                         #        yields well-differentiated modules; the turquoise module contains
                         #        the A. ricciae RVT homolog.
NETWORK_TYPE   <- "signed"
MIN_MODULE_SIZE<- 30
MERGE_CUTHEIGHT<- 0.25
DEEP_SPLIT     <- 2
MAX_BLOCK_SIZE <- 4000
N_THREADS      <- 10

MODULES_OF_INTEREST <- c("turquoise", "blue")   # <<< choose after inspecting the heatmap (section 9)

## ============================================================================
## 1. Libraries
## ============================================================================
# install.packages(c("tidyverse", "magrittr"))
# BiocManager::install(c("WGCNA", "DESeq2", "genefilter"))
library(tidyverse)
library(magrittr)
library(DESeq2)
library(WGCNA)

setwd(WORKDIR)
allowWGCNAThreads(nThreads = N_THREADS)

## ============================================================================
## 2. Load count matrix and define the sample design
## ============================================================================
data <- readr::read_delim(COUNTS_FILE, delim = "\t")
names(data)[1] <- "GeneId"
# Drop htseq-count summary rows if present (__no_feature, __ambiguous, ...)
data <- dplyr::filter(data, !grepl("^__", GeneId))

## Sample names: <ctrl|treatment>.<7|24><replicate letter>   e.g. ctrl.24b, treatment.7a
##   ctrl = uninfected control ; treatment = fungal infection
##   7 / 24 = hours post-infection (hpi) ; trailing letter = biological replicate
## Design: 4 groups (Control/Infected x 7h/24h), 3 replicates each = 12 samples.
sample_ids <- names(data)[-1]
meta_df <- data.frame(Sample = sample_ids) %>%
  mutate(
    Infection = ifelse(grepl("^ctrl", Sample), "Control", "Infected"),
    Time      = sub("^(ctrl|treatment)\\.(7|24)[a-z]$", "\\2", Sample),  # "7" or "24"
    Group     = paste0(Infection, "_", Time, "h")             # Control_7h, Infected_24h, ...
  ) %>%
  mutate(across(c(Infection, Time, Group), factor))
meta_df    # sanity-check the mapping before continuing

## ============================================================================
## 3. QC: distribution of counts per sample (spot outliers)
## ============================================================================
mdata <- data %>%
  tidyr::pivot_longer(cols = all_of(sample_ids)) %>%
  left_join(meta_df, by = c("name" = "Sample"))

ggplot(mdata, aes(x = name, y = value)) +
  geom_violin() +
  geom_point(alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Sample", y = "RNA-seq counts") +
  facet_grid(cols = vars(Group), drop = TRUE, scales = "free_x")

## ============================================================================
## 4. DESeq2 normalization
## ============================================================================
de_input <- as.matrix(data[, -1])
rownames(de_input) <- data$GeneId

dds <- DESeqDataSetFromMatrix(round(de_input), meta_df, design = ~ Group)
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts,
            file = paste0(OUT_PREFIX, "_normalized_counts.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

## ============================================================================
## 5. Build the matrix fed to WGCNA
## ----------------------------------------------------------------------------
## NOTE (important methodological choice):
##   The tutorial feeds a VARIANCE-STABILIZED, variance-filtered matrix to WGCNA
##   (recommended: WGCNA expects roughly homoscedastic, log-scale data).
##   USE_VST_FILTER = TRUE reproduces that. Setting it FALSE reproduces the
##   original run, which fed raw-scale DESeq-normalized counts instead.
## ============================================================================
if (USE_VST_FILTER) {
  wpn_vsd <- getVarianceStabilizedData(dds)            # VST (genes x samples)
  rv      <- genefilter::rowVars(wpn_vsd)
  keep    <- rv > quantile(rv, VAR_QUANTILE)           # top-variance genes
  expr_for_wgcna <- wpn_vsd[keep, ]
} else {
  expr_for_wgcna <- normalized_counts                  # original behaviour
}
message("Genes fed to WGCNA: ", nrow(expr_for_wgcna),
        " ; samples: ", ncol(expr_for_wgcna))

input_mat <- t(expr_for_wgcna)   # WGCNA wants samples x genes

## ============================================================================
## 6. Soft-thresholding power
## ============================================================================
powers <- c(1:10, seq(12, 40, by = 2))
sft <- pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2)); cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale-free topology fit, signed R^2",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"); abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
par(mfrow = c(1, 1))
# --> set PICKED_POWER in CONFIG from this plot, then continue.

## ============================================================================
## 7. Network construction and module detection
## ============================================================================
temp_cor <- cor
cor <- WGCNA::cor          # avoid namespace clash with stats::cor
netwk <- blockwiseModules(
  input_mat,
  power            = PICKED_POWER,
  networkType      = NETWORK_TYPE,
  deepSplit        = DEEP_SPLIT,
  pamRespectsDendro = FALSE,
  minModuleSize    = MIN_MODULE_SIZE,
  maxBlockSize     = MAX_BLOCK_SIZE,
  reassignThreshold = 0,
  mergeCutHeight   = MERGE_CUTHEIGHT,
  saveTOMs         = TRUE,
  saveTOMFileBase  = paste0(OUT_PREFIX, "_TOM"),
  numericLabels    = TRUE,
  verbose          = 3
)
cor <- temp_cor            # restore stats::cor

## ============================================================================
## 8. Module colors, dendrogram, and module table
## ============================================================================
mergedColors <- labels2colors(netwk$colors)

plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
table(mergedColors)

module_df <- data.frame(gene_id = names(netwk$colors),
                        colors  = mergedColors)
write_delim(module_df, file = paste0(OUT_PREFIX, "_gene_modules.txt"), delim = "\t")

## ============================================================================
## 9. Module eigengenes and module-sample heatmap
## ============================================================================
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order <- gsub("ME", "", names(MEs0))

MEs0$treatment <- rownames(MEs0)
mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(name = gsub("ME", "", name),
         name = factor(name, levels = module_order))

ggplot(mME, aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Module eigengenes across samples", y = "Modules", fill = "corr")

## ---- Module-trait correlation (infection status, timepoint, late infection) ----
traits <- meta_df %>%
  transmute(
    Sample,
    Infected     = as.integer(Infection == "Infected"),
    Time_h       = as.integer(as.character(Time)),
    Infected_24h = as.integer(Infection == "Infected" & as.character(Time) == "24")
  ) %>%
  column_to_rownames("Sample")
traits <- traits[rownames(input_mat), , drop = FALSE]     # align to expression rows

MEs <- orderMEs(moduleEigengenes(input_mat, mergedColors)$eigengenes)
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(input_mat))

textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits), yLabels = names(MEs), ySymbols = names(MEs),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
  setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1, 1),
  main = "Module-trait relationships"
)
par(mar = c(5, 4, 4, 2) + 0.1)

write.table(data.frame(module = rownames(moduleTraitCor), moduleTraitCor),
            file = paste0(OUT_PREFIX, "_module_trait_cor.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ============================================================================
## 10. Modules of interest: expression profiles + TOM edge list
## ============================================================================
rownames(module_df) <- module_df$gene_id
submod  <- subset(module_df, colors %in% MODULES_OF_INTEREST)
subexpr <- expr_for_wgcna[submod$gene_id, ]

data.frame(subexpr) %>%
  mutate(gene_id = rownames(.)) %>%
  pivot_longer(-gene_id) %>%
  mutate(module = module_df[gene_id, ]$colors) %>%
  ggplot(aes(x = name, y = value, group = gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module)) +
  labs(x = "sample", y = "normalized expression")

# Topological overlap for the selected genes -> Cytoscape/VisANT edge list
TOM <- TOMsimilarityFromExpr(t(subexpr), power = PICKED_POWER)
rownames(TOM) <- colnames(TOM) <- rownames(subexpr)

edge_list <- data.frame(TOM) %>%
  mutate(gene1 = rownames(.)) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(gene1 != gene2) %>%
  mutate(module1 = module_df[gene1, ]$colors,
         module2 = module_df[gene2, ]$colors)

write_delim(edge_list, file = paste0(OUT_PREFIX, "_edgelist.tsv"), delim = "\t")

################################################################################
## APPENDIX (optional): WGCNA on the DE-significant gene subset, step-by-step
##
## An alternative run used the significant-DEG subset (Ar_Sig_Exp_counts.txt)
## and the manual adjacency -> TOMdist -> cutreeDynamic route instead of
## blockwiseModules. Enable and adapt if you want this second analysis in the
## repo; parameters (power, minModuleSize) should be justified as above.
################################################################################
if (FALSE) {
  sig <- readr::read_delim("Ar_Sig_Exp_counts.txt", delim = "\t")
  datExpr <- t(data.frame(sig, row.names = 1))          # samples x genes
  gene.names <- colnames(datExpr)

  adj     <- adjacency(datExpr, power = PICKED_POWER, type = "signed")
  distTOM <- TOMdist(adj, TOMType = "signed")
  geneTree <- hclust(as.dist(distTOM), method = "average")

  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = distTOM,
                               deepSplit = 3, pamRespectsDendro = FALSE,
                               minClusterSize = 100)
  dynColors <- labels2colors(dynamicMods)
  plotDendroAndColors(geneTree, dynColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")

  # Write one gene list per module (drop unassigned 'grey')
  for (color in setdiff(unique(dynColors), "grey")) {
    module <- gene.names[dynColors == color]
    write.table(module,
                paste0(OUT_PREFIX, "_Sig_Exp_module_", color, ".txt"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
