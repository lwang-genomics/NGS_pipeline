#!/usr/bin/env python3
'''
Pipeline steps in RNA-seq
'''

import subprocess, sys
import utils

#===============
#===Tool Path===
#===============
FEATURECOUNTS = '/hpf/largeprojects/mdwilson/lib/subread-1.5.0-Linux-x86_64/bin/featureCounts'
HSAP_GTF='/hpf/largeprojects/mdwilson/genomes/hsap/gtf/gencode.v19.annotation.gtf'
MMUS_GTF='/hpf/largeprojects/mdwilson/genomes/mmus/gtf/gencode.vM11.annotation.gtf'
HSAP_SALMON="/hpf/largeprojects/mdwilson/genomes/hsap/salmon_index"
MMUS_SALMON="/hpf/largeprojects/mdwilson/genomes/mmus/salmon_index"
SALMON="/hpf/tools/centos6/salmon/0.13.1/bin/salmon"

#=====================
#===Pipeline Stages===
#=====================
@utils.formated_output
def star_se(star, gen_dir, prefix, threads, strand, star_genecounts=True):
    print('%s --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --runThreadN %d %s --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix %s --outWigType wiggle --outWigStrand %s'
                    %(star, gen_dir, prefix + '.trimmed.fastq.gz', threads, "--quantMode GeneCounts" if star_genecounts else ""
                      , prefix, strand) + '\n')

    subprocess.call('%s --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --runThreadN %d %s --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix %s --outWigType wiggle --outWigStrand %s'
                    %(star, gen_dir, prefix + '.trimmed.fastq.gz', threads, "--quantMode GeneCounts" if star_genecounts else "",
                      prefix, strand), shell=True)


@utils.formated_output
def star_pe(star, gen_dir, prefix1, prefix2, prefix, threads, strand, star_genecounts=True):
    print('%s --genomeDir %s --readFilesIn %s %s --readFilesCommand gunzip -c \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --runThreadN %d %s --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix %s --outWigType wiggle --outWigStrand %s'
                    %(star, gen_dir, prefix1 + '_paired.fastq.gz', prefix2 + '_paired.fastq.gz', threads,
                      "--quantMode GeneCounts" if star_genecounts else "", prefix, strand) + '\n')

    subprocess.call('%s --genomeDir %s --readFilesIn %s %s --readFilesCommand gunzip -c \
    --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --runThreadN %d %s --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix %s --outWigType wiggle --outWigStrand %s'
                    %(star, gen_dir, prefix1 + '_paired.fastq.gz', prefix2 + '_paired.fastq.gz', threads,
                      "--quantMode GeneCounts" if star_genecounts else "", prefix, strand), shell=True)


@utils.formated_output
def feature_counts(prefix, strand, threads, gtf_file, count_type="exon", pair_end=False):
    '''
    generate the counttable for a single library/sample
    '''
    strand_dict = {'forward': 1, 'reverse': 2, 'non': 0}
    print('%s %s -F GTF -t %s -T %d -s %d -g gene_name -a %s -o %s %s'
                    %(FEATURECOUNTS, '-p -B -C' if pair_end else '', count_type, threads, strand_dict[strand], gtf_file,
                      prefix + '_countable', prefix + '_roughuniq_sortedByName.bam' if pair_end else prefix +
                      '_roughuniq.bam') + '\n')

    subprocess.call('%s %s -F GTF -t %s -T %d -s %d -g gene_name -a %s -o %s %s'
                    %(FEATURECOUNTS, '-p -B -C' if pair_end else '', count_type, threads, strand_dict[strand], gtf_file,
                      prefix + '_countable', prefix + '_roughuniq_sortedByName.bam' if pair_end else prefix +
                      '_roughuniq.bam'), shell=True)


@utils.formated_output
def salmon(index, prefix1, prefix2, prefix, threads):
    print(f'{SALMON} quant -i {index} -l A -1 {prefix1}_paired.fastq.gz -2 {prefix2}_paired.fastq.gz '
                    f'-p {threads}  -o {prefix}_quant --validateMappings')
    subprocess.call(f'{SALMON} quant -i {index} -l A -1 {prefix1}_paired.fastq.gz -2 {prefix2}_paired.fastq.gz '
                    f'-p {threads} -o {prefix}_quant  --validateMappings', shell=True) ## note: it's highly recommended to use --validateMappings; but always show up the error message of unrecognized option.....



if __name__ == "__main__":
    print("Nothing to do, Please call from within rnaseq_stages.py")
    sys.exit()

