# NGS_pipeline: An All-in-One Modular Wrapper for RNA-seq, ChIP-seq, and ATAC-seq Workflows

**NGS_pipeline** is a streamlined, light-weight wrapper designed to simplify and automate the analysis of next-generation sequencing (NGS) data, including RNA-seq, ChIP-seq, and ATAC-seq. Built in Python, this tool provides flexible support for both single-end and paired-end reads, enabling users to perform trimming, alignment, quantification, peak calling, and quality control through intuitive commands. With minimal setup and clear structure, NGS_pipeline is ideal for both routine processing and reproducible research workflows.


## Installation

1. Clone the repository:

```bash
git clone https://github.com/lwang-genomics/NGS_pipeline.git
cd NGS_pipeline
```

2. Download default reference files (optional but recommended).  
Current default: hg38 & mm10

```bash
bash get_references.sh
```

3.  Install the pipeline locally:
```
pip install -e .
```
## External Dependencies

Please ensure the following tools are installed and accessible from your system PATH:

```text
## Common

Trimmomatic
SAMtools
Fastqc
MultiQC
wigToBigWig
Deeptools

## RNA-seq:
STAR
Salmon
Subread
Qualimap

## ChIP-seq & ATAC-seq:
BWA
MACS2
Aatqv

```

## I. RNA-Seq Processing Tool

This repository contains a modular and automated Python tool for processing RNA-seq data. It supports both single-end and paired-end reads with comprehensive support for trimming, alignment, quantification, and quality control. The tool offers flexibity to use both traditional mapping (STAR) and pseudo-alignment (Salmon).

### Features

- Automatic detection of read type (single-end or paired-end)  
- Trimming with **Trimmomatic**  
- Alignment using **STAR** or **Salmon**  
- Sorting and indexing of BAM files via **SAMtools**  
- Gene quantification with **featureCounts**   
- Quality control with **fastqc**, **Qualimap** (summarized by **MultiQC**) 
- Generation of strand-specific BigWig files  
- Detailed and reproducible logging at each processing stages 


### Example Usage

```bash
rna-seq SRR123456.R1.fq.gz SRR123456.R2.fq.gz\
  -sp hsap \
  -mq 10 \
  -st forward \
  -feature-level gene \
  --threads 4 &
```

### Options
```text
usage: rna-seq [-h] [-sp {hsap,mmus}] [-mq MAPQ] [-st {none,forward,reverse}] [-g GTF] [--pseudo] [--no-trim] [--feature-level {gene,exon}]
               [--keep-intermediate] [--genome-dir GENOME_DIR] [--threads THREADS]
               file1 [file2]

A light-weight RNA-seq preprocessing pipeline script. The tool now supports both single-end and paired-end RNA-seq and offers flexibility for both
traditional mapping (STAR) and psuedo alignment (Salmon)

positional arguments:
  file1                 Input FASTQ file (R1). For single-end reads, provide only this file.
  file2                 Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end.

options:
  -h, --help            show this help message and exit
  -sp {hsap,mmus}, --species {hsap,mmus}
                        Target species genome for alignment: 'hsap' (hg38) or 'mmus' (mm10). Default is 'hsap'.
  -mq MAPQ, --mapq MAPQ
                        Minimum MAPQ score to retain reads in the BAM file. Default is 5.
  -st {none,forward,reverse}, --strandness {none,forward,reverse}
                        Library strand-specific protocol. Choose from: none, forward, reverse. Default is reverse.
  -g GTF, --gtf GTF     Custom GTF file for feature counting. Overrides the default annotation.
  --pseudo              Use pseudo-alignment mode (e.g., Salmon) instead of alignment + featureCounts.
  --no-trim             Disable read trimming step. By default, trimming is performed.
  --feature-level {gene,exon}
                        Count features at the 'gene' or 'exon' level. Default is 'exon'.
  --keep-intermediate   Keep intermediate files (e.g., SAM, unfiltered BAM). Default is to remove them.
  --genome-dir GENOME_DIR
                        Optional: provide a custom genome directory path. Please structure the genome folder as defaut ones.
  --threads THREADS     Number of threads to use for multithreaded tools. Default is 4.
```
### Output
- sampleX*out.bam — Aligned and sorted BAM file
- sampleX_counts.txt — Gene count matrix
- sampleX*.bw — Strand-specific normalized BigWig signal
- sampleX_QC/ — Quality control reports


## II. ChIP-Seq Processing Tool

This pipeline provides an automated and lightweight solution for processing ChIP-seq data. It supports both single-end and paired-end sequencing reads and performs all major steps from raw read preprocessing to peak calling. The tool is optimized for simplicity and reproducibility.

## Features
- Supports both single-end and paired-end reads
- Trimming using **Trimmomatic**
- Alignment using **BWA**
- Conversion, filtering, and sorting of BAM files via **SAMtools**
- BigWig track generation with **deepTools** (bamCoverage)
- Peak calling using **MACS2**, with support for narrow and broad peak types
- Logging of each processing step with reproducible command outputs


### Example Usage

```bash
chip-seq SRR123456.R1.fq.gz\
  -sp hsap \
  -mq 10 \
  --peak-type narrow \
  --threads 4 &
```

### Options
```text
usage: chip-seq [-h] [-sp {hsap,mmus}] [-mq MAPQ] [--genome-dir GENOME_DIR] [--peak-type {narrow,broad}] [--no-trim] [--threads THREADS]
                [--keep-intermediate]
                file1 [file2]

A streamlined ChIP-seq preprocessing pipeline.

Features:
- Supports both single-end and paired-end reads.
- BWA-based alignment.
- Normalized bigWig signal track generation.
- Peak calling with MACS2 (narrow or broad peaks).

positional arguments:
  file1                 Input FASTQ file (R1). For single-end reads, provide only this file. Format: .(fastq|fq) or .(fastq|fq).gz
  file2                 Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end. Format: .(fastq|fq) or .(fastq|fq).gz

options:
  -h, --help            show this help message and exit
  -sp {hsap,mmus}, --species {hsap,mmus}
                        Target species genome: 'hsap' (hg38) or 'mmus' (mm10). Default is 'hsap'.
  -mq MAPQ, --mapq MAPQ
                        Minimum MAPQ score to retain reads in the BAM file. Default is 5.
  --genome-dir GENOME_DIR
                        Optional: provide a custom genome directory path. Please structure the genome folder as defaut ones.
  --peak-type {narrow,broad}
                        Type of peaks to call: narrow (default) or broad.
  --no-trim             Disable read trimming step. By default, trimming is performed.
  --threads THREADS     Number of threads. Default is 4.
  --keep-intermediate   Keep intermediate files.
```
### Output
- sampleX*.bam — Aligned and sorted BAM file
- sampleX*.bw — Library-size normalized BigWig signal
- sampleX_QC/ — Quality control reports
- sampleX\*Peak - peaks called by MACS2 



## III. ATAC-Seq Processing Tool

This pipeline offers a streamlined and reproducible solution for ATAC-seq data analysis. Designed for paired-end sequencing data, it performs essential steps from raw FASTQ files to peak calling and quality control. The pipeline prioritizes robustness, automation, and interpretability.

## Features

- Designed for paired-end reads
- Trimming using **Trimmomatic**
- Alignment using **BWA**
- Mitochondrial read filtering and MAPQ filtering using **SAMtools**
- BigWig track generation using **deepTools** (bamCoverage)
- Peak calling with **MACS2**, optimized for narrow peaks (with optional broad peak support)
- Quality control with **ATAQV** for ATAC-specific metrics
- Summary report generation using **MultiQC**
- Detailed logging for each stage to ensure reproducibility

### Example Usage

```bash
atac-seq SRR123456.R1.fq.gz SRR123456.R2.fq.gz\
  -sp hsap \
  -mq 10 \
  --threads 4 &
```

### Options
```text
usage: atac-seq [-h] [-sp {hsap,mmus}] [-mq MAPQ] [--genome-dir GENOME_DIR] [--peak-type {narrow,broad}] [--no-trim] [--threads THREADS]
                [--keep-intermediate]
                file1 file2

A streamlined ATAC-seq preprocessing pipeline.

Features:
- Supports paired-end reads.
- BWA-based alignment.
- Normalized bigWig signal track generation.
- Peak calling with MACS2 (default=narrowPeak).

positional arguments:
  file1                 Input FASTQ file (R1). Format: .(fastq|fq) or .(fastq|fq).gz
  file2                 Input FASTQ file (R2). Format: .(fastq|fq) or .(fastq|fq).gz

options:
  -h, --help            show this help message and exit
  -sp {hsap,mmus}, --species {hsap,mmus}
                        Target species genome: 'hsap' (hg38) or 'mmus' (mm10). Default is 'hsap'.
  -mq MAPQ, --mapq MAPQ
                        Minimum MAPQ score to retain reads in the BAM file. Default is 5.
  --genome-dir GENOME_DIR
                        Optional: provide a custom genome directory path. Please structure the genome folder as defaut ones.
  --peak-type {narrow,broad}
                        Type of peaks to call: narrow (default) or broad.
  --no-trim             Disable read trimming step. By default, trimming is performed.
  --threads THREADS     Number of threads. Default is 4.
  --keep-intermediate   Keep intermediate files.
```
### Output
- sampleX*.bam — Aligned and sorted BAM file
- sampleX*.bw — Library-size normalized BigWig signal
- sampleX\*Peak - peaks called by MACS2 
- sampleX_ataqv\* - ATAQV analysis
- sampleX_QC/ — Quality control reports





### Logging

Each major processing step is logged to a file with command and execution details (\*.log).

A test example:
```text
====================================================================================================
                                         [PIPELINE STARTED]
====================================================================================================
Start Time:          2025-06-01 08:27:08
Executed Command:    /opt/homebrew/Caskroom/mambaforge/base/bin/rna-seq SRR123456.R1.fq.gz SRR123456.R2.fq.gz -sp hsap -mq 10 -st forward --feature-level gene
====================================================================================================


====================================================================================================
                                 RNA-SEQ: PAIR-END MODE || CLASSIC
====================================================================================================


----------------------------------------------------------------------------------------------------
                                       STAGE: Quality Control
----------------------------------------------------------------------------------------------------
[COMMAND]       fastqc SRR123456.R1.fq.gz SRR123456.R2.fq.gz -t 4 -o SRR123456_QC

[STDOUT]
application/gzip
application/gzip
Analysis complete for SRR123456.R1.fq.gz
Analysis complete for SRR123456.R2.fq.gz


[STDERR]
Started analysis of SRR123456.R1.fq.gz
Approx 40% complete for SRR123456.R1.fq.gz
Approx 80% complete for SRR123456.R1.fq.gz
Started analysis of SRR123456.R2.fq.gz
Approx 40% complete for SRR123456.R2.fq.gz
Approx 80% complete for SRR123456.R2.fq.gz


[STATUS]        Stage 'Quality Control' completed successfully
----------------------------------------------------------------------------------------------------
.
.
.
.
.
.
----------------------------------------------------------------------------------------------------

====================================================================================================
                                        [PIPELINE COMPLETED]
====================================================================================================
End Time:            2025-06-01 08:58:19
====================================================================================================
```

### License

MIT License

### Acknowledgments

This pipeline integrates many excellent open-source bioinformatics tools. Credit goes to the developers of STAR, Salmon, BWA, SAMtools, Subread, Qualimap, UCSC tools and so on.



Please let me know if you have any questions or suggestions about my pipeline script!


