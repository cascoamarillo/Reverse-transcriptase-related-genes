# NcRVT-associated small RNA-seq - analysis pipeline

Analysis steps for the sequencing of small RNAs that co-purify with the native
NcRVT complex of *Neurospora crassa*, and their *in vitro* extension products.

## Overview

Native NcRVT complexes were purified from blasticidin S–induced *N. crassa*
(strain FGSC 2225) by sucrose density-gradient sedimentation. Co-purifying RNAs
were sequenced **before and after *in vitro* extension** by the NcRVT
template-independent polymerase, to identify NcRVT-associated / NcRVT-extendable
RNAs and to characterize the non-templated 3′ tails added by the enzyme.

**Samples (small RNA libraries, 150 bp paired-end, Illumina HiSeq 4000):**

| Sample | Fraction(s) | *In vitro* extension |
|--------|-------------|----------------------|
| `NcRVT_D4D5_noext` | D4+D5 (pooled peak) | none |
| `NcRVT_D4D5_NTP`   | D4+D5 (pooled peak) | ATP + CTP (NTPs) |
| `NcRVT_D4D5_dNTP`  | D4+D5 (pooled peak) | dATP + dCTP (dNTPs) |
| `NcRVT_D4`  | D4  (single fraction) | none |
| `NcRVT_D7`  | D7  (single fraction) | none |
| `NcRVT_D11` | D11 (single fraction) | none |

> Library prep: adapter-ligated with the NEBNext Multiplex Small RNA Library
> Prep Set for Illumina (requires 5′-phosphate / 3′-OH ends), size-selected for
> ~20–120 nt inserts on 10% PAGE, sequenced 150 bp paired-end (Genewiz-Azenta).

## Software

| Tool | Version | Purpose |
|------|---------|---------|
| Trimmomatic | 0.36 | adapter / low-quality trimming |
| FastQC | 0.11.4 | read QC |
| Bowtie | (Langmead 2009) | short-read alignment |
| SAMtools | — | BAM handling, size distribution |
| bioawk | — | read-length extraction |
| BEDTools | 2.29.2 | feature counting (`multicov`) |
| deepTools | 3.1.3 | coverage heatmaps (RPGC) |
| BWA-MEM | (Li & Durbin) | capture soft-clipped 3′ ends |
| SE-MEI | (dpryan79) | extract soft-clipped regions / breakpoints |
| qckitfastq | — | k-mer over-representation (`overrep_kmer`) |

**Reference genome:** *N. crassa* OR74A, FungiDB release 62 (https://fungidb.org).

## Pipeline

### 1. Adapter / quality trimming and QC

```bash
# Adapter and low-quality trimming
trimmomatic SE -phred33 sample_R1.fastq.gz sample_R1.trim.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:15

# Manual QC inspection
fastqc -o ./qc/ sample_R1.trim.fastq.gz
```

### 2. Alignment to the reference genome

```bash
bowtie -S -v 2 -a --best --strata \
    Nc_OR74A_release62_index \
    sample_R1.trim.fastq.gz sample.sam
samtools sort -o sample.sorted.bam sample.sam
samtools index sample.sorted.bam
```

### 3. Read-length (size) distribution

```bash
samtools view sample.sorted.bam \
  | bioawk -c sam '{ print length($seq) }' \
  | sort -n | uniq -c > sample.length_hist.txt
```

Preferential NcRVT extension is seen in the **40–120 nt** fraction; short RNAs
(<28 nt) correspond mostly to stress-induced tiRNAs/tRFs and qiRNAs.

### 4. Feature quantification

Counts across annotated *N. crassa* features (protein-coding CDS, rRNA, tRNA,
ncRNA), using the release-62 annotation:

```bash
bedtools multicov -bams sample.sorted.bam \
    -bed Nc_OR74A_release62_features.bed > sample.feature_counts.txt
```

### 5. Coverage heatmaps

```bash
# 1x depth (RPGC) normalized coverage
bamCoverage -b sample.sorted.bam --normalizeUsing RPGC \
    --effectiveGenomeSize <Nc_genome_size> -o sample.rpgc.bw

computeMatrix scale-regions -S *.rpgc.bw -R feature.bed -o matrix.gz
plotHeatmap -m matrix.gz -o heatmap.png
```

### 6. Non-templated 3′-tail analysis (soft-clipped reads)

The non-templated tails added by NcRVT appear as **soft-clipped** 3′ ends after
mapping. Reads were re-aligned with BWA-MEM to record the soft-clip coordinates,
and the clipped sequences were extracted with SE-MEI (minimum clip length 10 nt).

```bash
bwa mem Nc_OR74A_release62.fa sample_R1.trim.fastq.gz > sample.bwa.sam

# Extract soft-clipped regions, breakpoint positions, and full read sequences
samtools view -b sample.bwa.sam | samtools sort -o sample.bwa.bam
python SE-MEI/extractSoftclipped.py -l 10 sample.bwa.bam > sample.softclipped.fq.gz
```

### 7. k-mer over-representation of extension products

```bash
# observed/expected k-mer enrichment across soft-clipped sequences
overrep_kmer sample.softclipped.fq.gz -k 10 > sample.kmer_overrep.txt
```

The most over-represented 10-mers correspond to oligo(A) variants, consistent
with NcRVT adding A-rich non-templated tails.

## Data availability

- Raw + processed data: NCBI SRA, accession PRJNA1526284.
https://www.ncbi.nlm.nih.gov/sra/PRJNA1526284
