#!/usr/bin/env python3
'''
written in python 3.6.6

'''

import os,sys,time,argparse,subprocess,re,glob
from rnaseq_stages import *
from common_stages import *
import utils

# getting the current directory
this_dir = os.getcwd()


def main():
# Arguments:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file1',
                        help="input file name; can only take files from the current directory;\
                        accepted file format: sth.fastq, sth.fastq.gz, sth.sra")
    parser.add_argument('-f', '--file2', required=False,
                        help="only input file2, when -pm flag is used (i.e., pair-end)"
                             "an example name would be Something_r2.fastq.gz")
    parser.add_argument('-sp', '--species', required=False, default='hsap', choices=['hsap', 'mmus'],
                        help="choose the genome of species to align to")
    parser.add_argument('-t', '--threads', required=False, default=4, type=int,
                        help="number of threads to use")
    parser.add_argument('-st', '--strand', required=False, default='reverse', choices=['forward', 'reverse', 'non'],
                        help="specify strand-specific information")
    parser.add_argument('-mq', '--mapq', required=False, default=5, type=int,
                        help="mapq threshold used to filter out multiple mapping reads in bam files")
    parser.add_argument('-pm', '--pairend_mode', required=False, action="store_true",
                        help="If pair-end files, required to turn on this flag")
    parser.add_argument('-g', '--genome_dir', required=False,
                        help="the genome folder to use other than hg19 and mm10.(only support absolute paths)")
    parser.add_argument('-b', '--keep_raw_bam', required=False, action='store_true',
                        help="choose to keep the raw bam files")
    parser.add_argument('-c', '--count_type', required=False, default="exon",
                        help="specify the feature to count in featureCounts step; normal options are exon and gene")
    parser.add_argument('-sc', '--star_genecounts_off', required=False, action='store_true',
                        help="specify the feature to count in featureCounts step; normal options are exon and gene")
    parser.add_argument('-kw', '--keep_wiggle', required=False, action='store_true',
                        help="specify to keep the star outputted wiggle files")
    ####
    # psudoalignment flags:
    parser.add_argument('-ps', '--psuedo', required=False, action="store_true",
                        help="psuedo alignment using salmon")

    args = parser.parse_args()

# Check for valid arguments:
    if not os.path.exists(args.file1):
        print("Invalid File Name: File1 Does not Exist")
        sys.exit(1)

    if args.genome_dir:
        if not os.path.exists(args.genome_dir):
            print("Invalid Genome Folder: Folder Does not Exist")
            sys.exit(1)
        STAR = HSAP_STAR ## which star should be used for this??? maybe change to customized version
        GEN_DIR = args.genome_dir + '/star_index'
        CHROM_INFO = args.genome_dir + '/ChromInfo.txt'
        GTF = glob.glob(args.genome_dir + '/gtf/*.gtf')[0] # gtf is so messy!!! do not know which one to grab.... default: please leave only one .gtf file in the folder
        SALMON_INDEX = args.genome_dir + '/salmon_index'
    else:
        if args.species == 'hsap':
            STAR = HSAP_STAR
            GEN_DIR = HSAP_GEN_DIR
            CHROM_INFO = HSAP_CHROM_INFO
            GTF = HSAP_GTF
            SALMON_INDEX = HSAP_SALMON
        elif args.species == 'mmus':
            STAR = MMUS_STAR
            GEN_DIR = MMUS_GEN_DIR
            CHROM_INFO = MMUS_CHROM_INFO
            GTF = MMUS_GTF
            SALMON_INDEX = MMUS_SALMON

    # Check the starting format: // use fastq.gz as the standard starting format
    if re.search(r'\.sra$', args.file1):  ## note, re.match only matches from the beginning of the line; while re.search can match everywhere!
        sra_to_fastq_se(args.file1)
    elif re.search(r'\.fastq\.gz$', args.file1):
        pass
    elif re.search(r'\.fastq$', args.file1):
        subprocess.call('gzip %s' % args.file1, shell=True)
    else:
        print('ERROR: unsupported starting file format!!!')
        sys.exit(1)

# Check for single-end or paired-end:
    if not args.pairend_mode:
        # Run pipeline: (single end mode)
        print('='*30)
        print('RNA-SEQ: SINGLE-END MODE ')
        print('='*30)
        print('\n')

        ## start time:
        utils.print_time("Start")


        # get the prefix
        prefix = re.search(r'(.*?)\.', args.file1).group(1)   # use *? to do non-greedy match; and use back reference to get the first 1

        # start to run pipeline
        fastqc(prefix + '.fastq.gz', args.threads)
        trim_SE(prefix + '.fastq.gz', args.threads)
        star_se(STAR, GEN_DIR, prefix, args.threads, 'Unstranded' if (args.strand=='non') else 'Stranded',
                star_genecounts=not args.star_genecounts_off)
        os.renames(prefix + 'Aligned.sortedByCoord.out.bam', prefix + '.bam') # rename STAR_output files

        # index the raw bam files for qc
        index_bam(prefix + '.bam')

        roughuniq_bam(prefix + '.bam', args.threads, mapq=args.mapq)
        # sort_bam(prefix + '_roughuniq.bam', args.threads) # do not need to sort again...
        index_bam(prefix + '_roughuniq.bam')
        # may need to change the pattern if use different strand parameters
        # convert all the existing wig files to bw:
        wig_file = (re.match(prefix + r'.*Unique\.str[12]\.out\.wig', i).group() for i in os.listdir()
                    if re.match(prefix + r'.*Unique\.str[12]\.out\.wig', i))
        for i in wig_file: wiggleTobw(i, CHROM_INFO)
        feature_counts(prefix, args.strand, args.threads, GTF, count_type=args.count_type)

        # quality control for library complexity and bam files stats
        picard(prefix + '_roughuniq.bam')
        qualimap(prefix + '_roughuniq.bam', args.strand, GTF)

        multiqc()

        # generate the stats table: some from the multiqc folder, some by manual calculation
        generate_stats()
        print(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                        f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix}_stats_summary')
        print(f'awk \'$1=="{prefix}"{{print $6}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                        f'>>{prefix}_stats_summary')
        print(f'{SAMTOOLS} view -F 4 {prefix}.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                        f'>> {prefix}_stats_summary')
        print(f'{SAMTOOLS} view -q 60 {prefix}.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                        f'>> {prefix}_stats_summary')
        print(f'awk \'$1=="{prefix}"{{print $9}}\' ./multiqc_data/multiqc_featureCounts.txt | xargs echo -n -e "\\t"'
                        f'>>{prefix}_stats_summary')  # for exon mapping
        print(f'{SAMTOOLS} view {prefix}.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                        f'>> {prefix}_stats_summary')

        print(
            ['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                           f'\\tuniq_reads_perc\\texon_aligned_perc\\tmitocondria_reads_perc") <(awk \'BEGIN{{OFS="\\t"}}'
                           f'{{print $1,$2,$3,$4,$5/$2*100,$6/$2*100,$7/$2*100,$8/$2*100}}\' {prefix}_stats_summary) > '
                           f'{prefix}_stats_summary1'])
        print(f'mv {prefix}_stats_summary1 {prefix}_stats_summary')  # rename it
        print('\n')


        subprocess.call(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                        f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix}_stats_summary', shell=True)
        subprocess.call(f'awk \'$1=="{prefix}"{{print $6}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                        f'>>{prefix}_stats_summary', shell=True)
        subprocess.call(f'{SAMTOOLS} view -F 4 {prefix}.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                        f'>> {prefix}_stats_summary', shell=True)
        subprocess.call(f'{SAMTOOLS} view -q 60 {prefix}.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                        f'>> {prefix}_stats_summary', shell=True)
        subprocess.call(f'awk \'$1=="{prefix}"{{print $9}}\' ./multiqc_data/multiqc_featureCounts.txt | xargs echo -n -e "\\t"'
                        f'>>{prefix}_stats_summary', shell=True)  # for exon mapping
        subprocess.call(f'{SAMTOOLS} view {prefix}.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                        f'>> {prefix}_stats_summary', shell=True)

        subprocess.call(
            ['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                           f'\\tuniq_reads_perc\\texon_aligned_perc\\tmitocondria_reads_perc") <(awk \'BEGIN{{OFS="\\t"}}'
                           f'{{print $1,$2,$3,$4,$5/$2*100,$6/$2*100,$7/$2*100,$8/$2*100}}\' {prefix}_stats_summary) > '
                           f'{prefix}_stats_summary1'])
        # reorganize the table || note this can not run directly, because the default sh can not do parathesis;
        # so should explicitly specify using bash instead of sh... and then followed by the whole commands
        subprocess.call(f'mv {prefix}_stats_summary1 {prefix}_stats_summary', shell=True)  # rename it

        # clean up:
        subprocess.call(f'rm {prefix}.trimmed.fastq.gz {prefix}*.wig '
                        f'{prefix + ".bam*" if not args.keep_raw_bam else ""}', shell=True)

        print('='*30 + '\n')
        utils.print_time("Finish")
        print('COMPLETE SUCCESSFULLY!')
    else:
        # pair end mode
        if not args.file2:
            print("file2 for pairend pipeline is missing!")
            sys.exit(1)
        if not os.path.exists(args.file2):
            print("Invalid File Name: File2 Does not Exist")
            sys.exit(1)
        print('='*30)
        print(f'RNA-SEQ: PAIR-END MODE {"|| PSUEDO" if args.psuedo else ""}')
        print('='*30)
        print('\n')

        ## start time:
        utils.print_time("Start")

        # convert file2 also to fastq.gz format
        if re.search(r'\.sra$',
                     args.file2):  ## note, re.match only matches from the beginning of the line; while re.search can match everywhere!
            sra_to_fastq_se(args.file2)
        elif re.search(r'\.fastq\.gz$', args.file2):
            pass
        elif re.search(r'\.fastq$', args.file2):
            subprocess.call('gzip %s' % args.file2, shell=True)
        else:
            print('ERROR: unsupported starting file format!!!')
            sys.exit(1)

        prefix1 = re.search(r'(.*?)\.', args.file1).group(1)
        prefix2 = re.search(r'(.*?)\.', args.file2).group(1)
        # Need a better way to find the longest common subtring between two strings
        i = 0
        prefix = ''
        while prefix1[i] == prefix2[i]:
            if prefix1[i] == ".": ## make sure only take the common strings before the first dot
                break
            prefix += prefix1[i]
            i += 1
        # real start of pipeline
        fastqc(prefix1 + '.fastq.gz', args.threads)
        fastqc(prefix2 + '.fastq.gz', args.threads)
        trim_PE(prefix1, prefix2, args.threads)

        if args.psuedo :
            salmon(SALMON_INDEX, prefix1, prefix2, prefix, args.threads)
            # clean up:
            subprocess.call(f'rm {prefix}*paired.fastq.gz', shell=True)  # *_fastqc.zip
            
        else:
            star_pe(STAR, GEN_DIR, prefix1, prefix2, prefix, args.threads, 'Unstranded' if (args.strand == 'non') else 'Stranded',
                    star_genecounts=not args.star_genecounts_off)
            os.renames(prefix + 'Aligned.sortedByCoord.out.bam', prefix + '.bam') # rename STAR_output files
            # logics for post-processing pair-end bam files: filter out all single reads (mapq value > threshold),
            # filter out proper pairs, sort by read name, and then index

            #index raw bam files for qc
            index_bam(prefix + '.bam')

            roughuniq_bam(prefix + '.bam', args.threads, mapq=args.mapq)
            filter_pairs(prefix + '_roughuniq.bam', args.threads)
            os.renames(prefix + '_roughuniq_propaired.bam', prefix + '_roughuniq.bam')
            sort_bam(prefix + '_roughuniq.bam', args.threads, byname=True)
            os.renames(prefix + '_roughuniq_sorted.bam', prefix + '_roughuniq_sortedByName.bam')
            # sort_bam(prefix + '_roughuniq.bam', args.threads, byname=False) # do not need extra sorting
            os.rename(prefix + '_roughuniq.bam', prefix + '_roughuniq_sortedByPosition.bam')
            # turns out samtools index only accept position sorted bam files; so keep both version of bams for pair-end
            index_bam(prefix + '_roughuniq_sortedByPosition.bam')
            # another way is to use glob.glob (wildcards) to find matched files directly!!
            wig_file = (re.match(prefix + r'.*Unique\.str[12]\.out\.wig', i).group() for i in os.listdir()
                        if re.match(prefix + r'.*Unique\.str[12]\.out\.wig', i))
            for i in wig_file: wiggleTobw(i, CHROM_INFO)   ### should not do this!!! should only match the "prefix*wig" so do not have to do it when multiple files!!!!
            feature_counts(prefix, args.strand, args.threads, GTF, pair_end=True, count_type=args.count_type)

            # quality control for library complexity and bam files stats
            picard(prefix + '_roughuniq_sortedByPosition.bam')
            qualimap(prefix + '_roughuniq_sortedByName.bam', args.strand, GTF, pair_end=True)

            multiqc()

            # generate the stats table: some from the multiqc folder, some by manual calculation || using file1 to calculate..
            generate_stats()
            print(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix1}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                            f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix1}_stats_summary')
            print(
                f'awk \'$1=="{prefix1}"{{print 100-$3}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                f'>>{prefix1}_stats_summary')
            print(f'{SAMTOOLS} view -f 2 {prefix}.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                            f'>> {prefix1}_stats_summary')
            print(f'{SAMTOOLS} view -q 60 -f 2 {prefix}.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                            f'>> {prefix1}_stats_summary')
            print(f'awk \'$1=="{prefix1}"{{print $9}}\' ./multiqc_data/multiqc_featureCounts.txt | xargs echo -n -e "\\t"'
                            f'>>{prefix1}_stats_summary')  # for exon mapping
            print(
                f'{SAMTOOLS} view {prefix}.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                f'>> {prefix1}_stats_summary')

            print(
                ['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                               f'\\tuniq_reads_perc\\texon_aligned_perc\\tmitocondria_reads_perc") <(awk \'BEGIN{{OFS="\\t"}}'
                               f'{{print $1,$2,$3,$4,$5/($2*2)*100,$6/($2*2)*100,$7/($2*2)*100,$8/($2*2)*100}}\' {prefix1}_stats_summary) '
                               f'> {prefix1}_stats_summary1'])
            print(f'mv {prefix1}_stats_summary1 {prefix1}_stats_summary\n')  # rename it


            subprocess.call(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix1}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                            f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix1}_stats_summary', shell=True)
            subprocess.call(
                f'awk \'$1=="{prefix1}"{{print 100-$3}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                f'>>{prefix1}_stats_summary', shell=True)
            subprocess.call(f'{SAMTOOLS} view -f 2 {prefix}.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                            f'>> {prefix1}_stats_summary', shell=True)
            subprocess.call(f'{SAMTOOLS} view -q 60 -f 2 {prefix}.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                            f'>> {prefix1}_stats_summary', shell=True)
            subprocess.call(f'awk \'$1=="{prefix1}"{{print $9}}\' ./multiqc_data/multiqc_featureCounts.txt | xargs echo -n -e "\\t"'
                            f'>>{prefix1}_stats_summary', shell=True)  # for exon mapping
            subprocess.call(
                f'{SAMTOOLS} view {prefix}.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                f'>> {prefix1}_stats_summary', shell=True)

            subprocess.call(
                ['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                               f'\\tuniq_reads_perc\\texon_aligned_perc\\tmitocondria_reads_perc") <(awk \'BEGIN{{OFS="\\t"}}'
                               f'{{print $1,$2,$3,$4,$5/($2*2)*100,$6/($2*2)*100,$7/($2*2)*100,$8/($2*2)*100}}\' {prefix1}_stats_summary) '
                               f'> {prefix1}_stats_summary1'])
            # reorganize the table || note this can not run directly, because the default sh can not do parathesis;
            # so should explicitly specify using bash instead of sh... and then followed by the whole commands
            subprocess.call(f'mv {prefix1}_stats_summary1 {prefix1}_stats_summary', shell=True)  # rename it

            # clean up:
            subprocess.call(f'rm {prefix}*paired.fastq.gz {prefix + "*.wig " if not args.keep_wiggle else " "}'
                            f'{prefix + ".bam* " if not args.keep_raw_bam else ""}', shell=True)  #*_fastqc.zip

        print('='*30 + '\n')
        utils.print_time("Finish")
        print('COMPLETE SUCCESSFULLY!')


if __name__ == "__main__":
    main()
