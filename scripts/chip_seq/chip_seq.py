#!/usr/bin/env python3
"""
A streamlined ChIP-seq preprocessing pipeline.

Features:
- Supports both single-end and paired-end reads.
- BWA-based alignment.
- Normalized bigWig signal track generation.
- Peak calling with MACS2 (narrow or broad peaks).

"""

import argparse
import os
import subprocess
from pathlib import Path
import glob

from scripts.common import utils


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("file1", help="Input FASTQ file (R1). For single-end reads, provide only this file. Format: .(fastq|fq) or .(fastq|fq).gz")
    parser.add_argument("file2", nargs="?", help="Optional FASTQ file (R2) for paired-end reads. Leave blank for single-end. Format: .(fastq|fq) or .(fastq|fq).gz")
    parser.add_argument("-sp", "--species", choices=["hsap", "mmus"], default="hsap",
                        help="Target species genome: 'hsap' (hg38) or 'mmus' (mm10). Default is 'hsap'.")
    parser.add_argument('-mq', '--mapq', type=int, default=5,
                        help="Minimum MAPQ score to retain reads in the BAM file. Default is 5.")
    parser.add_argument("--genome-dir", type=Path, default=None,
                        help="Optional: provide a custom genome directory path. Please structure the genome folder as defaut ones.")
    parser.add_argument("--peak-type", choices=["narrow", "broad"], default="narrow",
                        help="Type of peaks to call: narrow (default) or broad.")
    parser.add_argument('--no-trim', action='store_true',
                        help="Disable read trimming step. By default, trimming is performed.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads. Default is 4.")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate files.")
    args = parser.parse_args()

    # Determine paired-end or single-end
    is_paired = args.file2 is not None

    if args.genome_dir:
        GENOME_DIR = args.genome_dir
        BWA_INDEX = GENOME_DIR / 'bwa_index'/ 'genome'
        CHROM_SIZES =list(GENOME_DIR.glob("*chrom.sizes"))[0]
    else:
        BASE_DIR = Path(__file__).resolve().parents[2]
        species_dir = 'hg38' if args.species == 'hsap' else 'mm10'
        GENOME_DIR = BASE_DIR / 'lib' / species_dir
        BWA_INDEX = GENOME_DIR / 'bwa_index'/ 'genome'
        CHROM_SIZES = GENOME_DIR / f'{species_dir}.chrom.sizes'

    # Paths and naming
    sample_name = utils.get_sample_name(args.file1)

    qc_dir = f"{sample_name}_QC"
    os.makedirs(qc_dir, exist_ok=True)

    log_file = utils.start_logging(sample_name)

    # Enhanced log header
    mode = "PAIR-END MODE" if is_paired else "SINGLE-END MODE"
    header = f"RNA-SEQ: {mode}"
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

    # Stage 3: Alignment with BWA
    sam_file = f"{sample_name}_aligned.sam"
    if is_paired:
        bwa_cmd = f"bwa mem -t {args.threads} {BWA_INDEX} {trimmed_R1} {trimmed_R2} > {sam_file}"
    else:
        bwa_cmd = f"bwa mem -t {args.threads} {BWA_INDEX} {trimmed_single} > {sam_file}"

    utils.log_stage("Alignment with BWA", bwa_cmd, log_file)

    # Stage 4: Convert, filter, sort and index BAM
    filtered_bam = f"{sample_name}_filtered_sorted.bam"
    filter_cmd = f"samtools view -bS -q {args.mapq} -@ {args.threads} {sam_file} | samtools sort -@ {args.threads} -o {filtered_bam} && samtools index {filtered_bam}"
    utils.log_stage("Convert, Filter and Sort BAM", filter_cmd, log_file)

    # Stage 5: Generate bigWig
    wig_cmd = f"bamCoverage -b {filtered_bam} -o {sample_name}.bw --binSize 10 --normalizeUsing CPM -p {args.threads}"
    utils.log_stage("Generate bigWig", wig_cmd, log_file)

    # Stage 6: Peak calling with MACS2
    if args.peak_type == "broad":
        peak_type_flag = "--broad --broad-cutoff 0.01"
    else:
        peak_type_flag = "-q 0.01"

    macs2_cmd = (
        f"macs2 callpeak -t {filtered_bam} -f {'BAMPE' if is_paired else 'BAM'} -n {sample_name} -g "
        f"{'hs' if args.species == 'hsap' else 'mm'} {peak_type_flag} "
        )
    utils.log_stage("MACS2 Peak Calling", macs2_cmd, log_file)

    # Stage 7: MultiQC summary
    multiqc_cmd = f"multiqc {qc_dir} -n {sample_name}_multiqc_report.html -o {qc_dir}"
    utils.log_stage("MultiQC Report", multiqc_cmd, log_file)

    # Cleanup
    if not args.keep_intermediate:
        files_to_remove = []
        if is_paired:
            files_to_remove = [trimmed_R1, trimmed_R2, unpaired_R1, unpaired_R2]
        else:
            files_to_remove = [trimmed_single]

        files_to_remove.append(sam_file)
        for f in files_to_remove:
            if f and os.path.exists(f):
                os.remove(f)

    utils.end_logging(log_file)

if __name__ == "__main__":
    main()



