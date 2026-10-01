# *N. crassa* mRNA-seq (wild-type vs *rvt*-null) - analysis pipeline

Analysis steps for poly(A) RNA-seq of blasticidin S–induced *Neurospora crassa*
wild-type and *rvt*-null strains, used to identify host genes co-regulated with
*rvt* (NcRVT, NCU09536).

## Overview

Four strains were profiled under blasticidin S induction and uninduced control
conditions, in biological replicates:

| Strain | Genotype | Background |
|--------|----------|-----------|
| FGSC 2489 | wild type (*mat A*) | reference |
| FGSC 4200 | wild type (*mat a*) | reference |
| H3 | *rvt*-null | derived from FGSC 2489 (*mat A*) |
| H2 | *rvt*-null | derived from FGSC 4200 (*mat a*) |

**Design:** 4 strains × {uninduced control, blasticidin S–induced} × 2–3
biological replicates = 20 libraries. Single-end, 50 bp, Illumina HiSeq 2000
(Brown University Genomics Core Facility).

> Library prep: poly(A)-selected RNA (MicroPoly(A)Purist, Ambion), strand-specific
> libraries with the Encore Complete RNA-Seq DR Multiplex System (NuGEN).

## Software

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.11.4 | read QC |
| cutadapt | 1.3 | 3′ adapter removal |
| FASTX-Toolkit (`fastq_quality_trimmer`) | — | quality trimming |
| TopHat (with Bowtie2) | — | spliced alignment |
| HTSeq (`htseq-count`) | — | per-gene read counts |
| DESeq (Bioconductor) | 1.x | differential expression (`nbinomTest`) |
| RColorBrewer / gplots / VennDiagram | — | heatmaps, PCA, Venn diagrams |

**Reference genome/annotation:** *N. crassa* OR74A v12 (FungiDB release 62):
`neurospora_crassa_or74a_12_supercontigs.fasta` +
`neurospora_crassa_or74a_12_transcripts.gtf` (TopHat) /
`neurospora_crassa_or74a_12_transcripts.gff3` (htseq-count).

## Pipeline

### 1. Read QC

```bash
fastqc -q -o ./raw/ -f fastq ~/illumina/Project_Ncrassa_lane4_140807/*.fastq.gz
```

### 2. Adapter removal (cutadapt)

```bash
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTCAGCATCTCGTATGCCGTCTTCTGCTTG

for f in ~/illumina/Project_Ncrassa_lane4_140807/*_R1_001.fastq.gz; do
    out=$(basename "${f%.fastq.gz}")_cutadapt.fastq.gz
    cutadapt -a "$ADAPTER" "$f" -o "$out"
done
```

*(Original single-file form:
`cutadapt -a $ADAPTER FER36_TTCAGC_L004_R1_001.fastq.gz -o FER36_TTCAGC_L004_R1_001_cutadapt.fastq.gz`)*

### 3. Quality trimming (FASTX-Toolkit)

Trim 3′ bases below Q20; discard reads shorter than 15 nt; Phred+33 encoding.

```bash
for f in *_cutadapt.fastq.gz; do
    zcat "$f" \
      | fastq_quality_trimmer -t 20 -l 15 -Q33 -z \
        -o "${f%_cutadapt.fastq.gz}_qtrim.fastq.gz"
done
```

| Flag | Meaning |
|------|---------|
| `-t 20` | trim 3′ nucleotides with quality < 20 |
| `-l 15` | minimum retained read length 15 nt |
| `-Q33` | Phred+33 quality encoding |
| `-z` | gzip output |

### 4. Spliced alignment (TopHat / Bowtie2)

```bash
GENOME=~/Ncgenome/neurospora_data/neurospora_crassa_or74a_12_supercontigs.fasta
GTF=~/Ncgenome/neurospora_data/neurospora_crassa_or74a_12_transcripts.gtf

tophat \
    -o ./tophat/Nc_2489_ctrl_rep1 \
    -p 12 \
    --library-type fr-firststrand \
    --no-convert-bam \
    --max-intron-length 400 \
    -G "$GTF" \
    "$GENOME" \
    Nc_2489_ctrl_rep1_qtrim.fastq.gz
```

Repeat per sample, changing `-o` and the input reads. Key parameters:

| Flag | Meaning |
|------|---------|
| `--library-type fr-firststrand` | strand-specific (NuGEN Encore) |
| `--max-intron-length 400` | *N. crassa* introns are short |
| `-G <gtf>` | gene-model-guided alignment |
| `--no-convert-bam` | keep SAM output |
| `-p 12` | threads |

### 5. Per-gene read counts (htseq-count)

Per-gene read counts were generated from the TopHat SAM alignments against the
OR74A gene models (`gene` features), one table per sample
(`Nc_<strain>_<cond>_rep<N>_geneID_count.txt`, columns: gene ID / raw count).

```bash
GFF3=~/Ncgenome/neurospora_data/neurospora_crassa_or74a_12_transcripts.gff3

htseq-count -q -s yes -t gene -i ID \
    ~/illumina/working/Nc_hiseq_140807/tophat/Nc_4200_ctrl_rep1.sam \
    "$GFF3" > Nc_4200_ctrl_rep1_geneID_count.txt
```

| Flag | Meaning |
|------|---------|
| `-s yes` | stranded counting (see note) |
| `-t gene` | feature type to count over |
| `-i ID` | GFF3 attribute used as the gene identifier |
| `-q` | suppress progress messages |

### 6. Combine per-sample counts into a matrix

For each comparison (e.g. 4200 vs H2, or 2489 vs H3), the per-sample count files
were sorted by gene ID and merged into one matrix used as DESeq input
(`Nc_4200_H2_geneID_counts.txt`, `Nc_2489_H3_geneID_count.txt`).

```bash
# Sort each count file by gene ID (keep header)
( head -1 "${in:=Nc_4200_ctrl_rep1_geneID_count.txt}" ; \
  tail -n +2 "$in" | sort -s -t , -k1,1 ) > my_counts_sorted.csv

# Merge sorted files on the gene-ID column (example, 3 files)
awk 'BEGIN{FS=","}
     FILENAME=="my_counts_sorted_1.csv" { dat[$1]=","$2 }
     FILENAME=="my_counts_sorted_2.csv" { print dat[$1]","$2 }' \
     my_counts_sorted_1.csv my_counts_sorted_2.csv \
  | paste my_counts_sorted_3.csv - > all_sorted.csv
```

### 7. Differential expression (DESeq, R)

Differential expression was computed with the **DESeq** package 

Condensed workflow (per combined matrix; example shown for 4200 vs H2):

```r
library("DESeq")

countTable <- read.table("Nc_4200_H2_geneID_counts.txt", header=TRUE)
rownames(countTable) <- countTable$geneID
countTable <- countTable[, -1]

# 3 ctrl + 3 bla per strain (2489/H3 matrices use 2 reps each)
condition <- factor(c("Nc4200_ctrl","Nc4200_ctrl","Nc4200_ctrl",
                      "Nc4200_bla","Nc4200_bla","Nc4200_bla",
                      "H2_ctrl","H2_ctrl","H2_ctrl",
                      "H2_bla","H2_bla","H2_bla"))

cds <- newCountDataSet(countTable, condition)
cds <- estimateSizeFactors(cds)
write.table(counts(cds, normalized=TRUE), "Nc_4200_H2_geneID_count_normalized.txt")
cds <- estimateDispersions(cds)

# Pairwise contrasts (induction within strain; genotype within condition)
res4200  <- nbinomTest(cds, "Nc4200_ctrl", "Nc4200_bla")
resH2    <- nbinomTest(cds, "H2_ctrl",     "H2_bla")
res_ctrl <- nbinomTest(cds, "Nc4200_ctrl", "H2_ctrl")
res_bla  <- nbinomTest(cds, "Nc4200_bla",  "H2_bla")

# Significant genes (Benjamini-Hochberg adjusted p < 0.05)
res4200Sig <- res4200[res4200$padj < 0.05, ]
write.csv(res4200Sig, "Nc4200Sig_padj_ctrl_vs_bla.csv")
```

Exploratory / visualization steps in the same script:
variance-stabilizing transformation (`varianceStabilizingTransformation`) →
sample-distance heatmap and PCA (`plotPCA`); top-gene expression heatmaps
(`heatmap.2`, gplots); and a 4-way Venn diagram of significant gene sets across
the four contrasts (`VennDiagram`).

## Data availability

- Raw + processed data: NCBI GEO, accession GSE348402 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348402
