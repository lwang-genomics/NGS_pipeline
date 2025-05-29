#!/usr/bin/env python3
'''

'''

import os,sys,time,argparse,subprocess,re,glob
from dnaseq_stages import *
from common_stages import *
import utils

def main():
# Arguments:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file1',
                        help="input r1 file name; can only take files from the current directory;\
                        accepted file format: sth.fastq, sth.fastq.gz, sth.sra")
    parser.add_argument('file2',
                        help="input r2 file name; can only take files from the current directory;\
                        accepted file format: sth.fastq, sth.fastq.gz, sth.sra")
    parser.add_argument('-sp', '--species', required=False, default='hsap', choices=['hsap', 'mmus'],
                        help="choose the genome of species to align to")
    parser.add_argument('-mq', '--mapq', required=False, default=5, type=int,
                        help="mapq threshold used to filter out multiple mapping reads in bam files")
    parser.add_argument('-qval', '--MACS2_qval', type=float, required=False, default=0.01,
                        help="q values when call peaks using macs2")
    parser.add_argument('-dp', '--dedup', required=False, action='store_true',
                        help="determine whether to dedup the duplicates from the library. Default is not to remove")
    parser.add_argument('-t', '--threads', required=False, default=4, type=int,
                        help="number of threads to use")
    parser.add_argument('-g', '--genome_dir', required=False,
                        help="the genome folder to use other than hg19 and mm10.(only support absolute paths)")
    parser.add_argument('-b', '--keep_raw_bam', required=False, action='store_true',
                        help="choose to keep the raw bam files")
    parser.add_argument('-rb', '--donot_remove_blacklist', required=False, action='store_true',
                        help="choose to not remove blacklist regions")

    args = parser.parse_args()

# Check for valid arguments:
    if not os.path.exists(args.file1):
        print("Invalid File1 Name: File1 Does not Exist")
        sys.exit(1)
    if not os.path.exists(args.file2):
        print("Invalid File2 Name: File2 Does not Exist")
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
    for i in [args.file1, args.file2]:
        if re.search(r'\.sra$', i):
            sra_to_fastq_se(i)
        elif re.search(r'\.fastq\.gz$', i):
            pass
        elif re.search(r'\.fastq$', i):
            subprocess.call('gzip %s' % i, shell=True)
        else:
            print('ERROR: unsupported starting file format!!!')
            sys.exit(1)

    # header
    print('=' * 30)
    print('ATAC-SEQ: PAIR-END')
    print('=' * 30)
    print('\n')

    ## start time:
    utils.print_time("Start")

    prefix1 = re.search(r'(.*?)\.', args.file1).group(1)
    prefix2 = re.search(r'(.*?)\.', args.file2).group(1)
    # find the prefix from prefix1 and prefix2
    i = 0
    prefix = ''
    while prefix1[i] == prefix2[i]:
        if prefix1[i] == ".":  ## make sure only take the common strings before the first dot
            break
        prefix += prefix1[i]
        i += 1

    # real start of pipeline
    fastqc(prefix1 + '.fastq.gz', args.threads)
    fastqc(prefix2 + '.fastq.gz', args.threads)
    trim_PE(prefix1, prefix2, args.threads, adapters=ATAC_ADAPTERS)
    bwa_pe(prefix1 + '_paired.fastq.gz', prefix2 + '_paired.fastq.gz', prefix, GENOME_FA, args.threads)
    sam_to_bam(prefix + '.sam', args.threads)

    #sort and index the raw bam files:
    sort_bam(prefix + '.bam', args.threads)
    index_bam(prefix + '_sorted.bam')

    roughuniq_bam(prefix + '_sorted.bam', args.threads, mapq=args.mapq)

    # remove black list regions
    if args.donot_remove_blacklist == False :
        remove_blacklist(prefix + '_sorted_roughuniq.bam', BLACKLIST)
        os.renames(prefix + '_sorted_roughuniq_blackout.bam', prefix + '_sorted_roughuniq.bam')

    filter_pairs(prefix + '_sorted_roughuniq.bam', args.threads)
    os.renames(prefix + '_sorted_roughuniq_propaired.bam', prefix + '_sorted_roughuniq.bam')

    # sort_bam(prefix + '_roughuniq.bam', args.threads, byname=False)
    os.renames(prefix + '_sorted_roughuniq.bam', prefix + '_roughuniq_sortedByPosition.bam')

    # turns out samtools index only accept position sorted bam files; so keep both version of bams for pair-end
    index_bam(prefix + '_roughuniq_sortedByPosition.bam')

    # remove mitocondria reads
    remove_mito_reads(prefix + '_roughuniq_sortedByPosition.bam')
    os.renames(prefix + '_roughuniq_sortedByPosition_filteredM.bam', prefix + '_roughuniq_sortedByPosition.bam')

    # resort the position sorted bam files to generate the name sorted bam files
    sort_bam(prefix + '_roughuniq_sortedByPosition.bam', args.threads, byname=True)
    os.renames(prefix + '_roughuniq_sortedByPosition_sorted.bam', prefix + '_roughuniq_sortedByName.bam')

    # whether to remove duplicates: seems picard can take both single and pair-end reads (sorted by position)
    if args.dedup:
        dedup(prefix + '_roughuniq_sortedByPosition.bam')

    bam_to_bw(prefix + '_roughuniq_sortedByPosition.bam', CHROM_INFO, args.threads, pe=True)

    # call peaks: for either QC or downstream analysis
    call_peak(prefix + '_roughuniq_sortedByName.bam', args.species,
              q_value=args.MACS2_qval, if_histone=True, pe=True)

    multiqc()

    #generate the stats table: some from the multiqc folder, some by manual calculation || using file1 to calculate..
    generate_stats()
    # print out the commands:
    print(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix1}"{{print $1"\\t"$14"\\t"100-$19}}\' '
          f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix1}_stats_summary')
    print(f'awk \'$1=="{prefix1}"{{print 100-$3}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
          f'>>{prefix1}_stats_summary')
    print(f'{SAMTOOLS} view -f 2 {prefix}_sorted.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                    f'>> {prefix1}_stats_summary')
    print(f'{SAMTOOLS} view -q 60 -f 2 {prefix}_sorted.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                    f'>> {prefix1}_stats_summary')
    print(f'{SAMTOOLS} view {prefix}_sorted.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                    f'>> {prefix1}_stats_summary')
    print(f'cat ${prefix}*broadPeak|wc -l|xargs echo -n -e "\\t" >> {prefix1}_stats_summary')
    print(['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                                   f'\\tuniq_reads_perc\\tmitocondria_reads_perc\\tpeak_num") <(awk \'BEGIN{{OFS="\\t"}}'
                                   f'{{print $1,$2,$3,$4,$5/($2*2)*100,$6/($2*2)*100,$7/($2*2)*100,$8}}\' {prefix1}_stats_summary) '
                                   f'> {prefix1}_stats_summary1'])
    print(f'mv {prefix1}_stats_summary1 {prefix1}_stats_summary')  # rename it
    print('\n')


    # execute the commands:
    subprocess.call(f'awk \'BEGIN{{FS="\\t"}}$1=="{prefix1}"{{print $1"\\t"$14"\\t"100-$19}}\' '
                    f'./multiqc_data/multiqc_fastqc.txt| xargs echo -n >> {prefix1}_stats_summary', shell=True)
    subprocess.call(f'awk \'$1=="{prefix1}"{{print 100-$3}}\' ./multiqc_data/multiqc_trimmomatic.txt | xargs echo -n -e "\\t"'
                    f'>>{prefix1}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view -f 2 {prefix}_sorted.bam -@ {args.threads} |wc -l|xargs echo -n -e "\\t" '
                    f'>> {prefix1}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view -q 60 -f 2 {prefix}_sorted.bam -@ {args.threads}|wc -l|xargs echo -n -e "\\t"'
                    f'>> {prefix1}_stats_summary', shell=True)
    subprocess.call(f'{SAMTOOLS} view {prefix}_sorted.bam -@ {args.threads} |grep chrM|wc -l| xargs echo -n -e "\\t" '
                    f'>> {prefix1}_stats_summary', shell=True)
    subprocess.call(f'cat ${prefix}*broadPeak|wc -l|xargs echo -n -e "\\t" >> {prefix1}_stats_summary', shell=True)

    subprocess.call(['bash', '-c', f'cat <(echo -e "sample\\traw_reads\\tduplicates_perc\\ttrimmed_perc\\taligned_perc'
                    f'\\tuniq_reads_perc\\tmitocondria_reads_perc\\tpeak_num") <(awk \'BEGIN{{OFS="\\t"}}'
                    f'{{print $1,$2,$3,$4,$5/($2*2)*100,$6/($2*2)*100,$7/($2*2)*100,$8}}\' {prefix1}_stats_summary) '
                    f'> {prefix1}_stats_summary1'])
                    # reorganize the table || note this can not run directly, because the default sh can not do parathesis;
                    # so should explicitly specify using bash instead of sh... and then followed by the whole commands
    subprocess.call(f'mv {prefix1}_stats_summary1 {prefix1}_stats_summary', shell=True)  # rename it

    # ataqv: use mkarv for the whole folder later...|| now ataqv only supports human and mouse (because have to input the species name)
    ataqv(prefix, args.species, args.threads)

    # clean up:
    subprocess.call(f'rm {prefix}*paired.fastq.gz {prefix}.bam '
                    f'{prefix + "_sorted.bam " + prefix + "_sorted.bam.bai" if not args.keep_raw_bam else ""} '
                    f'{prefix}.sam {prefix}*.bedGraph', shell=True)
    # if no *fastqc.zip, multiqc cannot capture them; choose to delete prefix_sorted.bam & prefix_sorted.bam.bai

    print('=' * 30 + '\n')
    utils.print_time("Finish")
    print('COMPLETE SUCCESSFULLY!')

if __name__ == "__main__":
    main()

