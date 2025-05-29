#!/usr/bin/env python3
"""
written in python 3.6.6

Refactored for clarity and simplified for shell-executable tools.
Now supports both single-end and paired-end RNA-seq samples.
"""

import os
import sys
import time
import argparse
import subprocess
import re
import glob

from ngs_pipeline.rnaseq_stages import *
from ngs_pipeline.common_stages import *
from ngs_pipeline.common import utils

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file1', help="Input FASTQ file (R1). For single-end reads, provide only this file.")
    parser.add_argument('file2', nargs='?', default=None,
                        help="Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end.")
    parser.add_argument('-sp', '--species', choices=['hsap', 'mmus'], default='hsap',
                        help="Target species genome for alignment: 'hsap' (human) or 'mmus' (mouse). Default is 'hsap'.")
    parser.add_argument('-mq', '--mapq', type=int, default=5,
                        help="Minimum MAPQ score to retain reads in the BAM file. Default is 5.")
    parser.add_argument('-st', '--strandness', choices=['none', 'forward', 'reverse'], default='reverse',
                        help="Library strand-specific protocol. Choose from: none, forward, reverse. Default is reverse.")
    parser.add_argument('-g', '--gtf', type=str, required=False,
                        help="Custom GTF file for feature counting. Overrides the default annotation.")
    parser.add_argument('--pseudo', required=False, action='store_true',
                        help="Use pseudo-alignment mode (e.g., Salmon) instead of alignment + featureCounts.")
    parser.add_argument('--no-trim', action='store_true',
                        help="Disable read trimming step. By default, trimming is performed.")
    parser.add_argument('--feature-level', choices=['gene', 'exon'], default='exon',
                        help="Count features at the 'gene' or 'exon' level. Default is 'exon'.")
    parser.add_argument('--keep-intermediate', action='store_true',
                        help="Keep intermediate files (e.g., SAM, unfiltered BAM). Default is to remove them.")
    parser.add_argument('--threads', type=int, default=4,
                        help="Number of threads to use for multithreaded tools. Default is 4.")
    args = parser.parse_args()

    BASE_DIR = Path(__file__).resolve().parents[1]
    species_dir = 'hg38' if args.species == 'hsap' else 'mm10'
    GENOME_DIR = BASE_DIR / 'lib' / species_dir
    STAR_INDEX = GENOME_DIR / 'star_index'
    SALMON_INDEX = GENOME_DIR / 'salmon_index'
    GTF_FILE = args.gtf if args.gtf else (GENOME_DIR / 'gtf' / ('gencode.v46.basic.annotation.gtf' if species_dir == 'hg38' else 'gencode.vM10.basic.annotation.gtf'))
    CHROM_SIZES = GENOME_DIR / f'{species_dir}.chrom.sizes'

    sample_name = utils.get_sample_name(args.file1)
    qc_dir = f"{sample_name}_QC"
    os.makedirs(qc_dir, exist_ok=True)

    log_file = utils.start_logging(sample_name)

    # Enhanced log header
    mode = "PAIR-END MODE" if args.file2 else "SINGLE-END MODE"
    protocol = "PSEUDO" if args.pseudo else "CLASSIC"
    header = f"RNA-SEQ: {mode} || {protocol}"
    log_file.write("\n" + "="*100 + "\n")
    log_file.write(f"{header.center(100)}\n")
    log_file.write("="*100 + "\n\n")
    log_file.flush()

    # Stage 1: Quality Control
    qc_cmd = f"fastqc {args.file1}{f' {args.file2}' if args.file2 else ''} -t {args.threads} -o {qc_dir}"
    utils.log_stage("Quality Control", qc_cmd, log_file)

    # Stage 2: Trimming (unless disabled)
    unpaired_R1 = unpaired_R2 = trimmed_R1 = trimmed_R2 = trimmed_single = None
    if args.no_trim:
        trimmed_files = (args.file1, args.file2) if args.file2 else (args.file1, None)
    else:
        if args.file2:
            trimmed_R1 = f"{sample_name}_R1_trimmed.fq.gz"
            trimmed_R2 = f"{sample_name}_R2_trimmed.fq.gz"
            unpaired_R1 = f"{sample_name}_R1_unpaired.fq.gz"
            unpaired_R2 = f"{sample_name}_R2_unpaired.fq.gz"
            trim_cmd = f"trimmomatic PE -threads {args.threads} {args.file1} {args.file2} {trimmed_R1} {unpaired_R1} {trimmed_R2} {unpaired_R2} -phred33 -gzip"
            trimmed_files = (trimmed_R1, trimmed_R2)
        else:
            trimmed_single = f"{sample_name}_trimmed.fq.gz"
            trim_cmd = f"java -jar trimmomatic-0.39.jar SE -threads {args.threads} {args.file1} {trimmed_single} -phred33 -gzip"
            trimmed_files = (trimmed_single, None)
        utils.log_stage("Trimming", trim_cmd, log_file)

    if args.pseudo:
        pseudo_cmd = f"salmon quant -i {SALMON_INDEX} -l A -r {trimmed_files[0]} {'-1 ' + trimmed_files[0] + ' -2 ' + trimmed_files[1] if args.file2 else ''} -p {args.threads} -o {sample_name}_salmon_output"
        utils.log_stage("Pseudo-alignment (Salmon)", pseudo_cmd, log_file)
    else:
        # Stage 3: Alignment with STAR
        aligned_prefix = f"{sample_name}_STAR"
        read_files_command = "--readFilesCommand zcat" if trimmed_files[0].endswith(".gz") else ""
        if args.file2:
            align_cmd = f"STAR --runThreadN {args.threads} --genomeDir {STAR_INDEX} --readFilesIn {trimmed_files[0]} {trimmed_files[1]} {read_files_command} --outFileNamePrefix {aligned_prefix}_ --outSAMtype BAM Unsorted --outWigType wiggle --outWigStrand Unstranded"
        else:
            align_cmd = f"STAR --runThreadN {args.threads} --genomeDir {STAR_INDEX} --readFilesIn {trimmed_files[0]} {read_files_command} --outFileNamePrefix {aligned_prefix}_ --outSAMtype BAM Unsorted --outWigType wiggle --outWigStrand Unstranded"
        utils.log_stage("Alignment (STAR)", align_cmd, log_file)

        # Stage 4: Filter, sort and index BAM
        unsorted_bam = f"{aligned_prefix}_Aligned.out.bam"
        sorted_bam = f"{sample_name}_filtered_sorted.bam"
        filter_sort_index_cmd = f"samtools view -b -q {args.mapq} {unsorted_bam} | samtools sort -@ {args.threads} -o {sorted_bam} && samtools index {sorted_bam}"
        utils.log_stage("Filter, Sort, and Index BAM", filter_sort_index_cmd, log_file)

        # Stage 5: Convert wiggle to bigWig
        wig_file = f"{aligned_prefix}_Signal.UniqueMultiple.str1.out.wig"
        bw_file = f"{sample_name}.bw"
        wig_to_bw_cmd = f"wigToBigWig {wig_file} {CHROM_SIZES} {bw_file}"
        utils.log_stage("Convert Wig to BigWig", wig_to_bw_cmd, log_file)

        # Stage 6: Feature counting
        count_output = f"{sample_name}_counts.txt"
        count_cmd = f"featureCounts -a {GTF_FILE} -o {count_output} -t {args.feature_level} -s {1 if args.strandness == 'forward' else 2 if args.strandness == 'reverse' else 0} -T {args.threads} {sorted_bam}"
        utils.log_stage("Feature Counting", count_cmd, log_file)

        # Stage 7: Picard QC
        picard_cmd = f"picard CollectRnaSeqMetrics I={sorted_bam} O={qc_dir}/{sample_name}_picard_metrics.txt REF_FLAT=ref_flat.txt STRAND={'SECOND_READ_TRANSCRIPTION_STRAND' if args.strandness == 'reverse' else 'FIRST_READ_TRANSCRIPTION_STRAND'} RIBOSOMAL_INTERVALS=rRNA.intervals"
        utils.log_stage("Picard RNA Metrics", picard_cmd, log_file)

        # Stage 8: Qualimap
        qualimap_cmd = f"qualimap rnaseq -outdir {qc_dir}/{sample_name}_qualimap -bam {sorted_bam} -gtf {GTF_FILE} -p {'strand-specific-reverse' if args.strandness == 'reverse' else 'strand-specific-forward' if args.strandness == 'forward' else 'non-strand-specific'} -outformat PDF"
        utils.log_stage("Qualimap RNA-seq QC", qualimap_cmd, log_file)

    # Stage 9: MultiQC summary
    multiqc_cmd = f"multiqc {qc_dir} -n {sample_name}_multiqc_report.html -o {qc_dir}"
    utils.log_stage("MultiQC Report", multiqc_cmd, log_file)

    # Cleanup intermediate files if requested
    if not args.keep_intermediate:
        for f in [trimmed_R1, trimmed_R2, unpaired_R1, unpaired_R2, trimmed_single, unsorted_bam, wig_file]:
            if f and os.path.exists(f):
                os.remove(f)

    utils.end_logging(log_file)

if __name__ == "__main__":
    main()