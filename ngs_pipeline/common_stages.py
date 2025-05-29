#!/usr/bin/env python3
'''
store all the common functions used for sequencing pipeline
'''
import subprocess, os
from ngs_pipeline.common import utils
#===============
#==Paths Info===
#===============
HSAP_STAR="/hpf/tools/centos6/star/2.5.1b/bin/STAR" #NOTE: mouse STAR and human STAR are not the same!!!!!!
MMUS_STAR="/hpf/largeprojects/mdwilson/external/STAR/bin/Linux_x86_64_static/STAR"
FASTQ_DUMP="/hpf/tools/centos6/sratoolkit/2.3.5-2/bin/fastq-dump"
FASTQC="/hpf/tools/centos6/fastqc/0.11.4/fastqc"
TRIMOMMATIC="/hpf/tools/centos6/trimmomatic/0.32/trimmomatic-0.32.jar" # trimmomatic parameters
ILLUMINACLIP="2:30:10"
WINDOWTRIM="4:15"
TRAILING='3'
LEADING='3'
MINLEN='36'  # have been changed to 20bp for LCL datasets, but the default is 36 bp; just keep it!!!
ADAPTERS="/hpf/tools/centos6/trimmomatic/source/Trimmomatic-0.32/adapters/TruSeq3-SE.fa"
ADAPTERS_PE="/hpf/tools/centos6/trimmomatic/source/Trimmomatic-0.32/adapters/TruSeq3-PE.fa"
ATAC_ADAPTERS="/hpf/tools/centos6/trimmomatic/source/Trimmomatic-0.32/adapters/NexteraPE-PE.fa"

#####################
SAMTOOLS="/hpf/tools/centos6/samtools/1.3.1/bin/samtools"
BEDTOOLS="/hpf/tools/centos6/bedtools/2.19.1/bin/bedtools"
MULTIQC="/hpf/largeprojects/mdwilson/external/Multiqc/multiqc"
JAVA="/hpf/tools/centos6/java/1.8.0_91/bin/java"
QUALIMAP="/hpf/largeprojects/mdwilson/external/qualimap_v2.2.1/qualimap"
FEATURECOUNTS="/hpf/largeprojects//mdwilson/external/subread-1.5.1-source/bin/featureCounts"
PICARD= "/hpf/tools/centos6/picard-tools/2.5.0/picard.jar"  # STAR parameters
HSAP_GEN_DIR="/hpf/largeprojects/mdwilson/genomes/hsap/star_index_2.5.1b"
MMUS_GEN_DIR="/hpf/largeprojects/mdwilson/genomes/mmus/star_index_2.4.1d_modified"
BAMCOVERAGE="/home/liangxi/rela/tools/deepTools_installed/bin/bamCoverage"
GENOMECOVERAGEBED="/hpf/tools/centos6/bedtools/2.27.1/bin/genomeCoverageBed"
HSAP_CHROM_INFO="/home/liangxi/mdwilson/genomes/hsap/ChromInfo_noRand.txt"
MMUS_CHROM_INFO="/hpf/largeprojects/mdwilson/genomes/mmus/ChromInfoNoScaffold.txt"
BEDGRAPHTOBIGWIG="/hpf/largeprojects/mdwilson/lib/bedGraphToBigWig"
wigToBigWig="/hpf/largeprojects/mdwilson/lib/wigToBigWig"


#==============
#===QC steps===
#==============

def fastqc(fastqgz, threads):
    '''
    quality control using fastqc
    '''
    print(' '.join([FASTQC, fastqgz, "-o", ".", "--noextract", "-t", str(threads)]) + '\n')
    subprocess.call([FASTQC, fastqgz, "-o", ".", "--noextract", "-t", str(threads)])



def picard(bam):
    print(' '.join([JAVA, '-jar', PICARD, 'EstimateLibraryComplexity', 'I=' + bam, 'O=' + bam.replace(".bam", ".txt")]) + '\n')
    subprocess.call([JAVA, '-jar', PICARD, 'EstimateLibraryComplexity', 'I=' + bam, 'O=' + bam.replace(".bam", ".txt")])



def qualimap(bam, strand, gtf, pair_end=False):
    '''
    quality control for bam files;
    input: bam file name, strand direction
    '''
    strand_dict = {"forward": "strand-specific-forward", 'reverse': 'strand-specific-reverse', 'non': 'non-strand-specific'}
    print('unset DISPLAY;%s rnaseq %s -bam %s -gtf %s -p %s --java-mem-size=8G -outdir %s'
                    %(QUALIMAP, '-pe -s' if pair_end else '', bam, gtf, strand_dict[strand], bam.replace(".bam", "")) + '\n')
    subprocess.call('unset DISPLAY;%s rnaseq %s -bam %s -gtf %s -p %s --java-mem-size=8G -outdir %s'
                    %(QUALIMAP, '-pe -s' if pair_end else '', bam, gtf, strand_dict[strand], bam.replace(".bam", "")), shell=True)



def multiqc():
    '''
    pipeline summary report by MultiQc
    '''
    ## NOTE!! have to use python 2.7.11 to run multiqc......
    python_version = '/hpf/tools/centos6/python/2.7.11/bin/python'
    print('module load python/2.7.6; %s %s . -f' %(python_version, MULTIQC) + '\n')
    subprocess.call('module load python/2.7.6; %s %s . -f' %(python_version, MULTIQC), shell=True)
    #seem that if call a new shell session in python 3, python3 path will be used by default; so need to module python2\
    #to run python2 scripts here.


#=====================
#===Pipeline Stages===
#=====================

def sra_to_fastq_se(sra):
    '''
    convert sra to fastq.gz format
    '''
    print(' '.join([FASTQ_DUMP, '--split-3', '--gzip', sra]) + '\n')
    subprocess.call([FASTQ_DUMP, '--split-3', '--gzip', sra])



def trim_SE(fastqgz, threads, adapters=ADAPTERS):
    '''
    Trim poor quility bases and/or adapter sequences from reads using Trimmomatic
    .fastq.gz -> trimmed.fastq.gz
    '''
    print(' '.join([JAVA, '-jar', '-Xmx8g', TRIMOMMATIC, 'SE', '-threads', str(threads), fastqgz, fastqgz.replace(".fastq.gz", ".trimmed.fastq.gz")
                    , 'ILLUMINACLIP:' + adapters + ':' + ILLUMINACLIP, 'TRAILING:'+TRAILING, 'LEADING:' + LEADING,
                     'MINLEN:' + MINLEN, 'SLIDINGWINDOW:' + WINDOWTRIM]) + '\n')
    subprocess.call([JAVA, '-jar', '-Xmx8g', TRIMOMMATIC, 'SE', '-threads', str(threads), fastqgz, fastqgz.replace(".fastq.gz", ".trimmed.fastq.gz")
                    , 'ILLUMINACLIP:' + adapters + ':' + ILLUMINACLIP, 'TRAILING:'+TRAILING, 'LEADING:' + LEADING,
                     'MINLEN:' + MINLEN, 'SLIDINGWINDOW:' + WINDOWTRIM])



def trim_PE(prefix1, prefix2, threads, adapters=ADAPTERS_PE):
    print('%s -jar -Xmx8g %s PE -phred33 -threads %s %s %s %s %s %s %s ILLUMINACLIP:%s:%s \
                    LEADING:%s TRAILING:%s SLIDINGWINDOW:%s MINLEN:%s' \
                    % (JAVA, TRIMOMMATIC, str(threads), prefix1 + '.fastq.gz', prefix2 + '.fastq.gz', prefix1 + '_paired.fastq.gz',
                      prefix1 + '_unpaired.fastq.gz', prefix2 + '_paired.fastq.gz', prefix2 + '_unpaired.fastq.gz', adapters,
                      ILLUMINACLIP, LEADING, TRAILING, WINDOWTRIM, MINLEN) + '\n')
    subprocess.call('%s -jar -Xmx8g %s PE -phred33 -threads %s %s %s %s %s %s %s ILLUMINACLIP:%s:%s \
                    LEADING:%s TRAILING:%s SLIDINGWINDOW:%s MINLEN:%s' \
                    % (JAVA, TRIMOMMATIC, str(threads), prefix1 + '.fastq.gz', prefix2 + '.fastq.gz', prefix1 + '_paired.fastq.gz',
                      prefix1 + '_unpaired.fastq.gz', prefix2 + '_paired.fastq.gz', prefix2 + '_unpaired.fastq.gz', adapters,
                      ILLUMINACLIP, LEADING, TRAILING, WINDOWTRIM, MINLEN), shell=True)


def sam_to_bam(sam, threads):
    print(f'{SAMTOOLS} view -S -b {sam} -@ {threads} > {sam.replace(".sam", ".bam")}' + '\n')
    subprocess.call(f'{SAMTOOLS} view -S -b {sam} -@ {threads} > {sam.replace(".sam", ".bam")}', shell=True)



def sort_bam(bam, threads, byname=False):
    print('%s sort %s --threads %s %s -o %s' %(SAMTOOLS, '-n' if byname else '', str(threads), bam,
                                                         bam.replace(".bam", "_sorted.bam")) + '\n')
    subprocess.call('%s sort %s --threads %s %s -o %s' %(SAMTOOLS, '-n' if byname else '', str(threads), bam,
                                                         bam.replace(".bam", "_sorted.bam")), shell=True)



def index_bam(bam):
    print(' '.join([SAMTOOLS, 'index', bam]) + '\n')
    subprocess.call([SAMTOOLS, 'index', bam])



def roughuniq_bam(bam, threads, mapq=5):
    #NOTE: '>' and '|' won't work in a list manner, have to use shell=True
    print('%s view -bq %s -@ %d %s > %s' %(SAMTOOLS, str(mapq), threads, bam, bam.replace(".bam", "_roughuniq.bam")) + '\n')
    subprocess.call('%s view -bq %s -@ %d %s > %s' %(SAMTOOLS, str(mapq), threads, bam, bam.replace(".bam", "_roughuniq.bam")),
                                               shell=True)



def filter_pairs(bam, threads):
    '''
    filter out the proper pairs in pair end bam files
    '''
    print('%s view -b -f 2 -@ %d %s > %s' %(SAMTOOLS, threads, bam, bam.replace(".bam", "_propaired.bam") + '\n'))
    subprocess.call('%s view -b -f 2 -@ %d %s > %s' %(SAMTOOLS, threads, bam, bam.replace(".bam", "_propaired.bam")), shell=True)



def bam_to_bw(bam, chrom_info, threads, pe=False):
    '''
    useful function to directly output normalized bigwig files;
    using RPM to normalize the file
    try if this looks good on browser [so macs2 do not normalize by rpm? correct??]
    [can also write a function to convert MACS outputed bedgraph to bw]
    improvements: keep a copy of raw bw file (without normalization)
    '''

    num = subprocess.check_output('%s view -@ %d %s | wc -l' %(SAMTOOLS, threads, bam), shell=True)
    if pe == True:
        num = int(num)/2

    # commands to generate raw bw files
    bam_to_bedGraph_cmd1 = '%s -bg -ibam %s  %s > %s' % (GENOMECOVERAGEBED, bam, '-pc' if pe else '', bam.replace('.bam', '.bedGraph'))

    bedGraph_to_bw_cmd1 = '%s %s %s %s' % (BEDGRAPHTOBIGWIG, bam.replace('.bam', '.bedGraph'), chrom_info, bam.replace('.bam', '_raw.bw'))

    print(bam_to_bedGraph_cmd1)
    print(bedGraph_to_bw_cmd1 + '\n')

    subprocess.call(bam_to_bedGraph_cmd1, shell=True)
    subprocess.call(bedGraph_to_bw_cmd1, shell=True)

    # commands to generate rpm normalized bw files
    bam_to_bedGraph_cmd2 = '%s -bg -ibam %s -scale %f %s > %s' % (GENOMECOVERAGEBED, bam, 1e6 / int(num),
                                                                 '-pc' if pe else '', bam.replace('.bam', '.bedGraph'))

    bedGraph_to_bw_cmd2 = '%s %s %s %s' % (BEDGRAPHTOBIGWIG, bam.replace('.bam', '.bedGraph'), chrom_info, bam.replace('.bam', '.bw'))

    print(bam_to_bedGraph_cmd2)
    print(bedGraph_to_bw_cmd2 + '\n')

    subprocess.call(bam_to_bedGraph_cmd2, shell=True)
    subprocess.call(bedGraph_to_bw_cmd2, shell=True)



def wiggleTobw(in_wig, chrom_sizes):
    print('%s %s %s %s' %(wigToBigWig, in_wig, chrom_sizes, in_wig.replace('.wig', '.bw')) + '\n')
    subprocess.call('%s %s %s %s' %(wigToBigWig, in_wig, chrom_sizes, in_wig.replace('.wig', '.bw')), shell=True)



def generate_stats():
    pass


