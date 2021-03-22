#!/usr/bin/env python3
'''

?? think about mapq values for BWA generated bam files

?? do we remove black regions from atac-seq datasets??

?? multiqc sometimes can sometimes can not recognize trimmomatic???

$$ errors with the cluster ataqv so have not implemented yet

NOTE: python 3.4 does not work; should use python 3.6.3 in cluster

ataqv also needs index....

$$ implement two default genome assembly: hg19 and mm10 which can be used for just type hsap/mmus ||
others should input the directory information

'''

import os,sys,time,argparse,subprocess,re,glob
from dnaseq_stages import *
from common_stages import *
import utils

# getting the current directory
this_dir = os.getcwd()


def main():
# Arguments:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file',
                        help="input file name; can only take files from the current directory;\
                            accepted file format: sth.fastq, sth.fastq.gz, sth.sra")
    parser.add_argument('-sp', '--species', required=False, default='hsap', choices=['hsap', 'mmus'],
                        help="choose the genome of species to align to")
    parser.add_argument('-in', '--input', required=False,
                        help="the input file (full name) used for peak calling")
    parser.add_argument('-his', '--histone_chip', required=False, action='store_true',
                        help="when dealing with histone chipseq, turn on this flag. Default is for TF chipseq")
    parser.add_argument('-mq', '--mapq', required=False, default=5, type=int,
                        help="mapq threshold used to filter out multiple mapping reads in bam files")
    parser.add_argument('-qval', '--MACS2_qval', type=float, required=False, default=0.01,
                        help="q values when call peaks using macs2")
    parser.add_argument('-atac', '--atac_seq', required=False, action='store_true',
                        help="turn on this flag if wants to run single end atac-seq using this script")
    parser.add_argument('-t', '--threads', required=False, default=4, type=int,
                        help="number of threads to use")
    parser.add_argument('-g', '--genome_dir', required=False,
                        help="the genome folder to use other than hg19 and mm10.(only support absolute paths)")
    parser.add_argument('-b', '--keep_raw_bam', required=False, action='store_true',
                        help="choose to keep the raw bam files")
    args = parser.parse_args()
    
# Check for valid arguments:
    if not os.path.exists(args.file):
        print("Invalid File Name: File Does not Exist")
        sys.exit(1)

    if args.input: # if input file is used
        if not os.path.exists(args.input):
            print("Invalid Input File Name: Input File Does not Exist")
            sys.exit(1)

    if args.MACS2_qval > 1 or args.MACS2_qval < 0:
        print("Invalid Q values")
        sys.exit(1)

    if args.genome_dir:
        if not os.path.exists(args.genome_dir):
            print("Invalid Genome Folder: Folder Does not Exist")
            sys.exit(1)
        GENOME_FA = glob.glob(args.genome_dir + '/bwa_index/*.fa')[0]
        BLACKLIST = args.genome_dir + '/blacklist.bed'
        CHROM_INFO = args.genome_dir + '/ChromInfo.txt'
    else:
        if args.species == 'hsap':
            GENOME_FA = HSAP_GENOME_FA
            BLACKLIST = HSAP_BLACKLIST
            CHROM_INFO = HSAP_CHROM_INFO
        elif args.species == 'mmus':
            GENOME_FA = MMUS_GENOME_FA
            BLACKLIST = MMUS_BLACKLIST
            CHROM_INFO = MMUS_CHROM_INFO

# Check the starting format: // use fastq.gz as the standard starting format
    if re.search(r'\.sra$', args.file):  ## note, re.match only matches from the beginning of the line; while re.search can match everywhere!
        sra_to_fastq_se(args.file)
    elif re.search(r'\.fastq\.gz$', args.file):
        pass
    elif re.search(r'\.fastq$', args.file):
        subprocess.call('gzip %s' % args.file, shell=True)
    else:
        print('ERROR: unsupported starting file format!!!')
        sys.exit(1)

    # header
    print('=' * 30)
    if not args.atac_seq:
        print('CHIP-SEQ: SINGLE-END')
    else:
        print('ATAC-SEQ: SINGLE-END')
    print('=' * 30)
    print('\n')

    ## start time:
    utils.print_time("Start")

    prefix = re.search(r'(.*?)\.', args.file).group(1)  # use *? to do non-greedy match; and use back reference to get the first 1

    # start to run pipeline
    fastqc(prefix + '.fastq.gz', args.threads)
    if args.atac_seq: # if atacseq, use another adapters to trim
        trim_SE(prefix + '.fastq.gz', args.threads, adapters=ATAC_ADAPTERS)
    else:
        trim_SE(prefix + '.fastq.gz', args.threads)
    bwa_se(prefix, GENOME_FA, args.threads)
    sam_to_bam(prefix + '.sam', args.threads)
    # sort and index the raw bam files:
    sort_bam(prefix + '.bam', args.threads)
    index_bam(prefix + '_sorted.bam')

    roughuniq_bam(prefix + '_sorted.bam', args.threads, mapq=args.mapq)
    # remove black list regions
    remove_blacklist(prefix + '_sorted_roughuniq.bam', BLACKLIST)
    index_bam(prefix + '_sorted_roughuniq_blackout.bam')

    if args.atac_seq:  # additional step for removing mitocondria reads for atac-seq
        remove_mito_reads(prefix + '_sorted_roughuniq_blackout.bam')
        os.renames(prefix + '_sorted_roughuniq_blackout_filteredM.bam', prefix + '_sorted_roughuniq_blackout.bam')

    bam_to_bw(prefix + '_sorted_roughuniq_blackout.bam', CHROM_INFO, args.threads)
    # de-duplicates: default is not to use this deduped files// just generate it here
    dedup(prefix + '_sorted_roughuniq_blackout.bam')

    # call peaks: for either QC or downstream analysis
    call_peak(prefix + '_sorted_roughuniq_blackout.bam', args.species,
              q_value=args.MACS2_qval, input=args.input, if_histone=args.histone_chip or args.atac_seq)

    multiqc()

    #generate the stats table: some from the multiqc folder, some by manual calculation
    generate_stats()
    print(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                    f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix}_stats_summary')
    print(f'awk \'$1=="{prefix}"{{print $6}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                    f'>>{prefix}_stats_summary')
    print(f'{SAMTOOLS} view -F 4 {prefix}_sorted.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                    f'>> {prefix}_stats_summary')
    print(f'{SAMTOOLS} view -q 60 {prefix}_sorted.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                    f'>> {prefix}_stats_summary')
    print(f'{SAMTOOLS} view {prefix}_sorted.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                    f'>> {prefix}_stats_summary')
    print(f'cat ${prefix}*broadPeak|wc -l|xargs echo -n -e "\\t" >> {prefix}_stats_summary')
    print(['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                                   f'\\tuniq_reads_perc\\tmitocondria_reads_perc\\tpeak_num") <(awk \'BEGIN{{OFS="\\t"}}'
                                   f'{{print $1,$2,$3,$4,$5/$2*100,$6/$2*100,$7/$2*100,$8}}\' {prefix}_stats_summary) > {prefix}_stats_summary1'])
    print(f'mv {prefix}_stats_summary1 {prefix}_stats_summary')  # rename it
    print('\n')



    subprocess.call(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                    f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix}_stats_summary', shell=True)
    subprocess.call(f'awk \'$1=="{prefix}"{{print $6}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                    f'>>{prefix}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view -F 4 {prefix}_sorted.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                    f'>> {prefix}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view -q 60 {prefix}_sorted.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                    f'>> {prefix}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view {prefix}_sorted.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                    f'>> {prefix}_stats_summary', shell=True)
    subprocess.call(f'cat ${prefix}*broadPeak|wc -l|xargs echo -n -e "\\t" >> {prefix}_stats_summary', shell=True)

    subprocess.call(['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                    f'\\tuniq_reads_perc\\tmitocondria_reads_perc\\tpeak_num") <(awk \'BEGIN{{OFS="\\t"}}'
                    f'{{print $1,$2,$3,$4,$5/$2*100,$6/$2*100,$7/$2*100,$8}}\' {prefix}_stats_summary) > {prefix}_stats_summary1'])
                    # reorganize the table || note this can not run directly, because the default sh can not do parathesis;
                    # so should explicitly specify using bash instead of sh... and then followed by the whole commands
    subprocess.call(f'mv {prefix}_stats_summary1 {prefix}_stats_summary', shell=True)  # rename it


    # clean up:
    subprocess.call(f'rm {prefix}.trimmed.fastq.gz {prefix}_sorted_roughuniq.bam '
                    f'{prefix}.sam {prefix}*.bedGraph {prefix}.bam {prefix + "_sorted.bam* " if not args.keep_raw_bam else ""}',
                    shell=True)
    # need to keep {prefix}_fastqc.zip, choose to delete prefix_sorted.bam, prefix_sorted.bam.bai

    print('=' * 30 + '\n')
    utils.print_time("Finish")
    print('COMPLETE SUCCESSFULLY!')


if __name__ == '__main__':
    main()

