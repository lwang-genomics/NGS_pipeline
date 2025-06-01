# NGS_pipeline: RNA-Seq Processing Tool

This repository contains a modular and automated Python tool for processing RNA-Seq data. It supports both single-end and paired-end reads and covers trimming, alignment, quantification, and quality control.

## Features

- **Automatic detection** of read type (single-end or paired-end)  
- **Trimming** with Trimmomatic  
- **Alignment** using STAR  
- **Sorting and indexing** of BAM files via SAMtools  
- **Gene quantification** with featureCounts  
- **RNA-Seq metrics** from Picard  
- **Qualimap** summary statistics  
- **Optional output of strand-specific BigWig files**  
- **Logging** of commands and stages for reproducibility  

## File Structure
NGS_pipeline/
├── scripts/
│   ├── rna_seq.py        # Main executable script
│   └── utils.py          # Logging and helper functions
├── lib/                  # Genome indices, annotation files, etc.
└── README.md

## Installation

1. Clone the repository:

```bash
git clone https://github.com/yourusername/NGS_pipeline.git
cd NGS_pipeline
```
2.	(Optional) Install dependencies in a virtual environment:
```bash
pip install -r requirements.txt
```

3.	Build and install the package locally:

## Usage
```bash
rna-seq \
  --sample-name SRR123456 \
  --input R1.fq.gz R2.fq.gz \
  --genome-dir lib/hg38/star_index \
  --gtf lib/hg38/annotation.gtf \
  --ref-flat lib/hg38/refFlat.txt \
  --rrna-intervals lib/hg38/rRNA.intervals \
  --chrom-sizes lib/hg38/chrom.sizes \
  --threads 4
```

## Options
	•	--sample-name : Sample identifier
	•	--input : Input FASTQ files (1 or 2)
	•	--genome-dir : STAR genome index directory
	•	--gtf : GTF annotation file
	•	--ref-flat : RefFlat file for Picard
	•	--rrna-intervals : Interval file for rRNA
	•	--chrom-sizes : Chromosome sizes file for BigWig conversion
	•	--strandness : none, forward, or reverse (default: reverse)
	•	--keep-intermediate : Keep intermediate files (default: False)

## Output
	•	Trimmed FASTQ files
	•	STAR-aligned BAM file (sorted + indexed)
	•	Gene count matrix
	•	Quality control reports (Picard, Qualimap)
	•	Optional BigWig signal tracks

## Logging

Each major processing step is logged to a file with timestamps and command details.

## License

MIT License

## Author

Liangxi Wang – [lwang@ipmc.cnrs.fr]

Let me know if you have any questions!



