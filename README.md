# RNA-seq

This repository contains a complete RNA-seq analysis pipeline designed for SLURM environments.

## Contents
- Pre-processing
- Differential Gene Expresion Analysis

## Pre-procesing Step-by-Step Description

### How to Run

1. Modify the variables inside `pipeline.sh` with your paths and sample names.
2. Submit to SLURM:
   ```bash
   sbatch RNAseq_PIPELINE.sh
   ```

### 1. Quality Control (FastQC)
Runs a quality control check on the raw FASTQ files.

Command:

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o output_dir

```
Arguments:

- sample.fastq.gz: Input FASTQ file
- -o: Output directory

### 2. Alignment (STAR)

Aligns paired-end reads to the reference genome using STAR. Outputs include:
- Sorted BAM file
- Gene count matrix

Command:

```bash
STAR \
  --genomeDir STAR_index_dir \
  --sjdbGTFfile annotation.gtf \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 4 \
  --outFileNamePrefix output_prefix \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts

```
- --genomeDir:	Path to STAR genome index (e.g. built for hg38)
- --sjdbGTFfile:	GTF annotation file used for guided alignment
- --readFilesIn:	Input paired-end FASTQ files
- --readFilesCommand:	Decompression method (e.g., zcat for .gz)
- --runThreadN:	Number of threads for parallel processing
- --outFileNamePrefix:	Prefix for all output files
- --outSAMtype: BAM SortedByCoordinate	Output sorted BAM file by coordinate
- --outTmpDir:	Temporary directory for STAR to use
- --outWigStrand: Unstranded	Type of strand specificity in wiggle files (Unstranded, Forward, or Reverse)
- --quantMode: GeneCounts	Outputs read counts per gene


### 3. Index BAM File (SAMtools)
Indexes the sorted BAM file for downstream tools.

Command: 

```bash
samtools index Aligned.sortedByCoord.out.bam

```


### 4. Generate Signal Tracks (deepTools)
Creates a normalized bigWig file for visualization.

Command:

```bash
bamCoverage \
  --bam input.bam \
  --outFileName output.bw \
  --effectiveGenomeSize 2913022398 \
  --outFileFormat bigwig \
  --binSize 1 \
  --normalizeUsing RPGC
```

Arguments:

- --bam: Input BAM file (cleaned, filtered, and deduplicated)
- --outFileName: Output bigWig file name
- --effectiveGenomeSize: Effective genome size (e.g. 2913022398 for human hg38)
- --outFileFormat: Output format, typically bigwig
- --binSize: Bin size for signal aggregation (e.g. 1 bp resolution)
- --normalizeUsing: Normalization method (e.g. RPGC = Reads Per Genomic Content)

## Differential Gene Expresion Analysis
