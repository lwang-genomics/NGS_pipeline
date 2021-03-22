#!/usr/bin/env python3
'''
Pipeline stages specific to DNA-seq
'''

import subprocess, sys
import utils

#===============
#==Paths Info===
#===============
BWA="/hpf/tools/centos6/bwa/0.7.8/bin/bwa"
HSAP_GENOME_FA="/hpf/largeprojects/mdwilson/genomes/hsap/bwa_index/Homo_sapiens.GRCh37.68.dna.chromosomes.fa"
MMUS_GENOME_FA="/hpf/largeprojects/mdwilson/genomes/mmus/bwa_index/Mus_musculus.GRCm38.68.dna.chromosomes.fa"
MACS2="/hpf/tools/centos6/python/2.7.11/bin/macs2"
JAVA="/hpf/tools/centos6/java/1.8.0_91/bin/java"
PICARD= "/hpf/tools/centos6/picard-tools/2.5.0/picard.jar"

HSAP_BLACKLIST='/hpf/largeprojects/mdwilson/genomes/hsap/blacklist.bed'
MMUS_BLACKLIST='/hpf/largeprojects/mdwilson/genomes/mmus/mm10.blacklist.bed'
BEDTOOLS='/hpf/tools/centos6/bedtools/2.27.1/bin/bedtools'
SAMTOOLS="/hpf/tools/centos6/samtools/1.3.1/bin/samtools"


#=====================
#===Pipeline Stages===
#=====================
@utils.formated_output
def bwa_se(prefix, genome_fa, threads):
    '''
    DNA-seq aligner
    '''
    print(f'{BWA} mem -M -t {threads} {genome_fa} {prefix}.trimmed.fastq.gz > {prefix}.sam' + '\n')
    subprocess.call(f'{BWA} mem -M -t {threads} {genome_fa} {prefix}.trimmed.fastq.gz > {prefix}.sam', shell=True)


@utils.formated_output
def bwa_pe(fastqgz1, fastqgz2, prefix, genome_fa, threads):
    '''
    DNA-seq aligner for pair-end
    '''
    cmd = f'{BWA} mem -M -t {threads} {genome_fa} {fastqgz1} {fastqgz2} > {prefix}.sam'
    print(cmd)
    subprocess.call(cmd, shell=True)


@utils.formated_output
def call_peak(bam, species, q_value=0.01, input=None, if_histone=False, pe=False):
    '''
    default mode is without input bam file, and this will be used for quality control purpose
    (check if the peak number is reasonable)
    very default parameters, only work on TF chipseq; for histone chip, should use --nomodel
    '''

    species_dict = {"hsap": 'hs', "mmus": 'mm'}
    print(  f'module load python/2.7.6; {MACS2} callpeak -t {bam} -f {"BAMPE" if pe else "BAM"} {"-c " + input if input else ""} '
            f'-g {species_dict[species]} -n {bam.replace(".bam","")} '
            f'{"--nomodel --broad --broad-cutoff " + str(q_value) if if_histone else "-q " + str(q_value)}' + '\n')

    subprocess.call(f'module load python/2.7.6; {MACS2} callpeak -t {bam} -f {"BAMPE" if pe else "BAM"} '
                    f'{"-c " + input if input else ""} -g {species_dict[species]} -n {bam.replace(".bam","")} '
                    f'{"--nomodel --broad --broad-cutoff " + str(q_value) if if_histone else "-q " + str(q_value)}'
    ,shell=True)


@utils.formated_output
def dedup(bamfile):
    print(f'{JAVA} -jar -Xmx8g {PICARD} MarkDuplicates INPUT={bamfile} REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE={bamfile.replace(".bam","_dedup.metrics")} \
    CREATE_INDEX=true OUTPUT={bamfile.replace(".bam", "_dedup.bam")}' + '\n')

    subprocess.call(f'{JAVA} -jar -Xmx8g {PICARD} MarkDuplicates INPUT={bamfile} REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE={bamfile.replace(".bam","_dedup.metrics")} \
    CREATE_INDEX=true OUTPUT={bamfile.replace(".bam", "_dedup.bam")}', shell=True)


@utils.formated_output
def remove_blacklist(bam, blacklist_region):
    '''
    remove reads from bam files
    '''
    print(f'{BEDTOOLS} intersect -abam {bam} -b {blacklist_region} -v > '
                    f'{bam.replace(".bam","_blackout.bam")}' + '\n')

    subprocess.call(f'{BEDTOOLS} intersect -abam {bam} -b {blacklist_region} -v > '
                    f'{bam.replace(".bam","_blackout.bam")}', shell=True)


@utils.formated_output
def remove_mito_reads(bamfile):
    print(f'{SAMTOOLS} idxstats {bamfile} | cut -f 1 | grep -v MT | xargs {SAMTOOLS} view -b {bamfile} > '
                    f'{bamfile.replace(".bam","_filteredM.bam")}')
    subprocess.call(f'{SAMTOOLS} idxstats {bamfile} | cut -f 1 | grep -v MT | xargs {SAMTOOLS} view -b {bamfile} > '
                    f'{bamfile.replace(".bam","_filteredM.bam")}', shell=True)


@utils.formated_output
def ataqv(prefix, species, threads):
    species_dict = {'hsap':'human', 'mmus':'mouse'}
    tss_dict = {'hsap':'/hpf/largeprojects/mdwilson/liangxi/projects/rela/public_data/genomes/hg19/gencode.v19.gene.tss.annotation.bed',
                'mmus':'/hpf/largeprojects/mdwilson/liangxi/projects/rela/public_data/genomes/mm10/gencode.vM11.gene.annotation.tss.bed'}
    cmd = f'module load ataqv/1.0.0; ataqv --threads {threads} {species_dict[species]} ' \
          f'{prefix}_sorted.bam --peak-file {prefix}*broadPeak --tss-file {tss_dict[species]}'
    print(cmd)
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    print("Nothing to do, Please call from within dna-seq python scripts")
    sys.exit()
