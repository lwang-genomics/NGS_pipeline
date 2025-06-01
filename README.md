# NGS_pipeline: RNA-Seq Processing Tool

This repository contains a modular and automated Python tool for processing RNA-Seq data. It supports both single-end and paired-end reads with comprehensive support for trimming, alignment, quantification, and quality control. The tool offers flexibity to use both traditional mapping (STAR) and pseudo-alignment (Salmon).

## Features

- Automatic detection of read type (single-end or paired-end)  
- Trimming with Trimmomatic  
- Alignment using **STAR** or **Salmon**  
- Sorting and indexing of BAM files via SAMtools  
- Gene quantification with **featureCounts**   
- Quality control with **fastqc**, **Qualimap** (summarized by **MultiQC**) 
- Generation of strand-specific **BigWig** files  
- Detailed and reproducible logging at each processing stages 

## Installation

1. Clone the repository:

```bash
git clone https://github.com/lwang-genomics/NGS_pipeline.git
cd NGS_pipeline
```

2.	Install the pipeline locally:
```
pip install .
```
## External Dependencies

Please ensure the following tools are installed and accessible from your system PATH:

STAR
Salmon
Trimmomatic
samtools
subread
qualimap
wigToBigWig


## Example Usage

```bash
rna-seq SRR123456.R1.fq.gz SRR123456.R2.fq.gz\
  -sp hsap \
  -mq 10 \
  -st forward \
  -feature-level gene \
  --threads 4 &
```

## Options
usage: rna-seq [-h] [-sp {hsap,mmus}] [-mq MAPQ] [-st {none,forward,reverse}] [-g GTF] [--pseudo] [--no-trim] [--feature-level {gene,exon}]
               [--keep-intermediate] [--threads THREADS]
               file1 [file2]

A light-weight RNA-seq preprocessing pipeline script. Now supports both single-end and paired-end RNA-seq samples Supports both classic (STAR)
alignment and psuedo alignment (salmon)

positional arguments:
  file1                 Input FASTQ file (R1). For single-end reads, provide only this file.
  file2                 Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end.

options:
  -h, --help            show this help message and exit
  -sp {hsap,mmus}, --species {hsap,mmus}
                        Target species genome for alignment: 'hsap' (human) or 'mmus' (mouse). Default is 'hsap'.
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
  --threads THREADS     Number of threads to use for multithreaded tools. Default is 4.

## Output
	•	sampleX*out.bam — Aligned and sorted BAM file
	•	sampleX_counts.txt — Gene count matrix
	•	sampleX*.bw — Strand-specific normalized BigWig signal
	•	sampleX_QC/ — Quality control reports

## Logging

Each major processing step is logged to a file with command and execution details (\*.log).


## License

MIT License

## Acknowledgments

This pipeline integrates many excellent open-source bioinformatics tools. Credit goes to the developers of STAR, Salmon, Trimmomatic, SAMtools, Subread, Qualimap, UCSC tools and so on.

Let me know if you have any questions or suggestions about my pipeline script!




