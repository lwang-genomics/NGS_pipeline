GENOME_DIR = Path("lib/hg38")
#!/usr/bin/env python3
"""
A light-weight RNA-seq preprocessing pipeline script.
Now supports both single-end and paired-end RNA-seq samples
Supports both classic (STAR) alignment and psuedo alignment (salmon)

"""

import os
import sys
import time
import argparse
import subprocess
import re
import glob
from pathlib import Path 

from scripts.rnaseq_stages import *
from scripts.common_stages import *
from scripts.common import utils

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genome-dir', type=Path, default=args.genome_dir,
                        help='Custom genome directory path. Must match expected folder structure.')
    parser.add_argument('file1', help="Input FASTQ file (R1). For single-end reads, provide only this file.")
    parser.add_argument('file2', nargs='?', default=None,
                        help="Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end.")
    parser.add_argument('-sp', '--species', choices=['hsap', 'mmus'], default='hsap',
                        help="Target species genome for alignment: 'hsap' (hg38) or 'mmus' (mm10). Default is 'hsap'.")
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

    BASE_DIR = Path(__file__).resolve().parents[2]
    species_dir = 'hg38' if args.species == 'hsap' else 'mm10'
    args.genome_dir = BASE_DIR / 'lib' / species_dir
    STAR_INDEX = args.genome_dir / 'star_index'
    SALMON_INDEX = args.genome_dir / 'salmon_index'
    GTF_FILE = args.gtf if args.gtf else (args.genome_dir / 'gtf' / ('gencode.v46.basic.annotation.gtf' if species_dir == 'hg38' else 'gencode.vM10.basic.annotation.gtf'))
    CHROM_SIZES = args.genome_dir / f'{species_dir}.chrom.sizes'

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
            trim_cmd = f"trimmomatic PE -threads {args.threads} {args.file1} {args.file2} {trimmed_R1} {unpaired_R1} {trimmed_R2} {unpaired_R2} -phred33 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25"
            trimmed_files = (trimmed_R1, trimmed_R2)
        else:
            trimmed_single = f"{sample_name}_trimmed.fq.gz"
            trim_cmd = f"trimmomatic SE -threads {args.threads} {args.file1} {trimmed_single} -phred33 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25"
            trimmed_files = (trimmed_single, None)
        utils.log_stage("Trimming", trim_cmd, log_file)

    if args.pseudo:
        pseudo_cmd = f"salmon quant -i {SALMON_INDEX} -l A -r {trimmed_files[0]} {'-1 ' + trimmed_files[0] + ' -2 ' + trimmed_files[1] if args.file2 else ''} -p {args.threads} -o {sample_name}_salmon_output"
        utils.log_stage("Pseudo-alignment (Salmon)", pseudo_cmd, log_file)
    else:
        # Stage 3: Alignment with STAR
        aligned_prefix = f"{sample_name}_STAR"
        read_files_command = "--readFilesCommand gzcat" if trimmed_files[0].endswith(".gz") else ""
        if args.file2:
            align_cmd = f"STAR --runThreadN {args.threads} --genomeDir {STAR_INDEX} --readFilesIn {trimmed_files[0]} {trimmed_files[1]} {read_files_command} --outFileNamePrefix {aligned_prefix}_ --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM"
        else:
            align_cmd = f"STAR --runThreadN {args.threads} --genomeDir {STAR_INDEX} --readFilesIn {trimmed_files[0]} {read_files_command} --outFileNamePrefix {aligned_prefix}_ --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM"
        utils.log_stage("Alignment (STAR)", align_cmd, log_file)

        # Stage 4: Filter, sort and index BAM
        unsorted_bam = f"{aligned_prefix}_Aligned.sortedByCoord.out.bam"
        sorted_bam = f"{sample_name}_filtered_sorted.bam"
        filter_sort_index_cmd = f"samtools view -b -q {args.mapq} {unsorted_bam} | samtools sort -@ {args.threads} -o {sorted_bam} && samtools index {sorted_bam}"
        utils.log_stage("Filter, Sort, and Index BAM", filter_sort_index_cmd, log_file)

        # Stage 5: Convert wiggle to bigWig
        wig_cmds = []
        for strand in ['str1', 'str2']:
            wig_file = f"{aligned_prefix}_Signal.Unique.{strand}.out.wig"
            bw_file = f"{sample_name}.{strand}.bw"
            wig_cmds.append(f"wigToBigWig {wig_file} {CHROM_SIZES} {bw_file}")

        combined_cmd = " && ".join(wig_cmds)
        utils.log_stage("Convert Wig to BigWig (str1 & str2)", combined_cmd, log_file)

        # Stage 6: Feature counting
        bam_file = f"{sample_name}_filtered_sorted.bam"
        counts_file =f"{sample_name}_counts.txt"
        paired_flag = "-p" if utils.is_paired_end_bam(bam_file) else ""
        if args.strandness == "none":
            strand_option = "-s 0"
        elif args.strandness == "forward":
            strand_option = "-s 1"
        else:  # "reverse"
            strand_option = "-s 2"

        featurecounts_cmd = (
            f"featureCounts {paired_flag} -T {args.threads} "
            f"{strand_option} -a {GTF_FILE} -o {counts_file} {bam_file}"
        )
        utils.log_stage("Feature Counting", featurecounts_cmd, log_file)

        # Stage 7: Qualimap
        strand_option = {
            "reverse": "strand-specific-reverse",
            "forward": "strand-specific-forward",
            "none": "non-strand-specific"
        }[args.strandness]

        qualimap_cmd = (
            f"qualimap rnaseq "
            f"-outdir {qc_dir}/{sample_name}_qualimap "
            f"-bam {sorted_bam} "
            f"-gtf {GTF_FILE} "
            f"-p {strand_option} "
            f"--java-mem-size=4G "
            f"-outformat PDF"
        )
        utils.log_stage("Qualimap RNA-seq QC", qualimap_cmd, log_file)

    # Stage 8: MultiQC summary
    multiqc_cmd = f"multiqc {qc_dir} -n {sample_name}_multiqc_report.html -o {qc_dir}"
    utils.log_stage("MultiQC Report", multiqc_cmd, log_file)

    # Cleanup intermediate files if requested
    if not args.keep_intermediate:
        files_to_remove = [trimmed_R1, trimmed_R2, unpaired_R1, unpaired_R2, trimmed_single, unsorted_bam]
        files_to_remove += glob.glob(f"{aligned_prefix}_Signal.Unique.*.out.wig")

        for f in files_to_remove:
            if f and os.path.exists(f):
                os.remove(f)

    utils.end_logging(log_file)

if __name__ == "__main__":
    main()