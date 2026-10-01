# rvt - reverse transcriptase-related genes at the intersection of stress response pathways

Analysis code, pipelines, and processing notes for the manuscript:

> **Reverse transcriptase-related genes at the intersection of stress response pathways**
> Yushenova I.A., Rodriguez F., Wilkinson M.E., Metzger B.M., Gladyshev E.A., Zhang F., Arkhipova I.R.

This repository documents the transcriptomic and co-expression analyses behind
the paper. *rvt* genes are a distinct, ancient class of cellular reverse
transcriptases that are strongly induced by environmental and translational
stress; the analyses here characterize the host genes and pathways co-regulated
with *rvt* in the fungus *Neurospora crassa* and the bdelloid rotifer
*Adineta ricciae*.

## Contents

| Analysis | Files |
|----------|-------|
| Small RNA-seq of NcRVT-associated RNAs and their in vitro extension products (*N. crassa*) | `sRNA-seq_analysis.md` |
| Poly(A) mRNA-seq of blasticidin-induced *N. crassa* wild-type vs *rvt*-null strains (DESeq differential expression) | `mRNA-seq_analysis.md`, `deseq_commands.R` |
| WGCNA co-expression network of the *A. ricciae* fungal-infection time course | `README_WGCNA.md`, `wgcna_Ar.R` |

Each analysis has its own documented script / notes with the exact commands,
tool versions, and parameters.

## Analyses at a glance

**Neurospora small RNA-seq** — RNAs co-purifying with the native NcRVT complex
(blasticidin S–induced *N. crassa* FGSC 2225), sequenced before and after in
vitro extension with dNTPs/NTPs to identify NcRVT-associated and
NcRVT-extendable RNAs and the non-templated 3′ tails added by the enzyme.
Pipeline: Trimmomatic → Bowtie → feature quantification (bedtools) → soft-clip
3′-tail analysis (BWA-MEM / SE-MEI / qckitfastq).

**Neurospora mRNA-seq** — poly(A) RNA-seq of two wild-type strains (FGSC 2489,
FGSC 4200) and their isogenic *rvt*-null mutants (H3, H2), under blasticidin S
induction vs uninduced control, in biological replicates. Pipeline:
cutadapt → FASTX quality trimming → TopHat/Bowtie2 → htseq-count → DESeq.

**Adineta ricciae WGCNA** — re-analysis of the fungal-infection RNA-seq time
course from Nowell et al. 2024, adding weighted gene co-expression network
analysis (WGCNA). A signed network (power 20) resolves co-expression modules;
the turquoise module contains the *A. ricciae* RVT homolog.

## Data availability

- Neurospora small RNA-seq: NCBI SRA, PRJNA152628
- Neurospora mRNA-seq: NCBI GEO, GSE348402 
- Adineta ricciae RNA-seq: from Nowell et al. 2024, *Nat. Commun.* 15:5787,
  https://doi.org/10.1038/s41467-024-49919-1 (BioProject RJEB39927 https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB39927)

## Software

R (DESeq / DESeq2, WGCNA, tidyverse) and standard command-line RNA-seq tools
(Trimmomatic, cutadapt, FASTX-Toolkit, FastQC, Bowtie/Bowtie2, TopHat, HTSeq,
SAMtools, BEDTools, deepTools, BWA, ViennaRNA). Exact versions and commands are
in each subfolder's documentation.

## Citation

If you use this code, please cite the manuscript above. Method references:

- Langfelder P, Horvath S (2008). WGCNA. *BMC Bioinformatics* 9:559. https://doi.org/10.1186/1471-2105-9-559
- Anders S, Huber W (2010). Differential expression analysis for sequence count data (DESeq). *Genome Biol.* 11:R106.
- Love MI, Huber W, Anders S (2014). DESeq2. *Genome Biol.* 15:550.
- Nowell RW et al. (2024). Bdelloid rotifers deploy horizontally acquired biosynthetic genes against a fungal pathogen. *Nat. Commun.* 15:5787.
