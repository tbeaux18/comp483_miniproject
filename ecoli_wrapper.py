#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Timothy Baker
@version: 1.0.0

ecoli_wrapper.py

This module is a python wrapper that pulls the FASTA files from specific FTP URLs from NCBI.


6. The assembled genome in RefSeq for E. coli K-12 (NC_000913) has 4140 CDS and 89 tRNAs annotated.
Write to the log file the discrepancy (if any) found. For instance, if my Prokka annotation predicted 4315 CDS and 88 tRNA’s,
I would write, Prokka found 175 additional CDS and 1 less tRNA than the RefSeq in assembly HM27.
[You’ll likewise write out the number of bp for the other 3 strains.]

7. Now that we know where the genes are located, we can see how these genes are transcribed.
    Use TopHat & Cufflinks to map the reads of a specific strain to the genome of the strain and quantify their expression, respectively.

9. Using Cuffnorm, normalize all 4 of your transcriptomes.
    Parse the output of Cuffnorm (the Simple-table gene attributes format file) such that you create a sorted file,
    with the highest expressed gene first, for each strain and write this to file, e.g. HM27_normalized_sorted.tsv.

Dependencies:
    python3 v > 3.5
    Prokka
    tophat2
        tophat
        bowtie2
    sratools
        prefetch
        fastq-dump
    cufflinks
        cufflinks
        cuffmerge
        cuffnorm
    Biopython

"""

import os
import sys
import argparse
import subprocess
import logging
from Bio import SeqIO

 #  _____     _______ _    _  _____
 # |  __ \ /\|__   __| |  | |/ ____|
 # | |__) /  \  | |  | |__| | (___
 # |  ___/ /\ \ | |  |  __  |\___ \
 # | |  / ____ \| |  | |  | |____) |
 # |_| /_/    \_\_|  |_|  |_|_____/
 #

# (FASTA FTP LINK, FEATURE COUNT FTP LINK)
HM27_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz')

HM46_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz')

HM65_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz')

HM69_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz')



FIRST_LAST_PATH = './Timothy_Baker'
if os.path.isdir('./Timothy_Baker'):
    os.chdir(FIRST_LAST_PATH)
else:
    os.mkdir(FIRST_LAST_PATH)
    os.chdir(FIRST_LAST_PATH)

CURRENT_DIR = os.getcwd()
print(CURRENT_DIR)
sys.exit()
FASTA_DIR_PATH = './ncbi_fasta'

HM27_FASTA = './ncbi_fasta/HM27_FASTA.fna'
HM46_FASTA = './ncbi_fasta/HM46_FASTA.fna'
HM65_FASTA = './ncbi_fasta/HM65_FASTA.fna'
HM69_FASTA = './ncbi_fasta/HM69_FASTA.fna'

HM27_GFF_FILE = './hm27_index.gff'
HM46_GFF_FILE = './hm46_index.gff'
HM65_GFF_FILE = './hm65_index.gff'
HM69_GFF_FILE = './hm69_index.gff'

HM27_BAM = './hm27_tophat/accepted_hits.bam'
HM46_BAM = './hm46_tophat/accepted_hits.bam'
HM65_BAM = './hm65_tophat/accepted_hits.bam'
HM69_BAM = './hm69_tophat/accepted_hits.bam'

HM27_SORTED_BAM = './hm27_tophat/accepted_hits.sorted.bam'
HM46_SORTED_BAM = './hm46_tophat/accepted_hits.sorted.bam'
HM65_SORTED_BAM = './hm65_tophat/accepted_hits.sorted.bam'
HM69_SORTED_BAM = './hm69_tophat/accepted_hits.sorted.bam'

HM27_FASTQ_1 = './hm27_sra/SRR1278956_1.fastq'
HM27_FASTQ_2 = './hm27_sra/SRR1278956_2.fastq'
HM46_FASTQ_1 = './hm46_sra/SRR1278960_1.fastq'
HM46_FASTQ_2 = './hm46_sra/SRR1278960_2.fastq'
HM65_FASTQ_1 = './hm65_sra/SRR1283106_1.fastq'
HM65_FASTQ_2 = './hm65_sra/SRR1283106_2.fastq'
HM69_FASTQ_1 = './hm69_sra/SRR1278963_1.fastq'
HM69_FASTQ_2 = './hm69_sra/SRR1278963_2.fastq'

MERGED_GTF = './merged_ecoli/merged.gtf'

 #  _      ____   _____  _____ ______ _____
 # | |    / __ \ / ____|/ ____|  ____|  __ \
 # | |   | |  | | |  __| |  __| |__  | |__) |
 # | |   | |  | | | |_ | | |_ |  __| |  _  /
 # | |___| |__| | |__| | |__| | |____| | \ \
 # |______\____/ \_____|\_____|______|_|  \_\
 #

LOGGER = logging.getLogger(__name__)

LOGGER.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')

FILE_HANDLER = logging.FileHandler("working.log")

FILE_HANDLER.setFormatter(FORMATTER)

LOGGER.addHandler(FILE_HANDLER)



def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-t', '--threads', help='Integer declaring how many threads to run')

    return parser.parse_args()



 #  ______       _____ _______
 # |  ____/\    / ____|__   __|/\
 # | |__ /  \  | (___    | |  /  \
 # |  __/ /\ \  \___ \   | | / /\ \
 # | | / ____ \ ____) |  | |/ ____ \
 # |_|/_/    \_\_____/   |_/_/    \_\
 #

def wget_gunzip_fasta(ftp_list):
    """ grabs the zipped fasta and feature txt files for the specific strains
        listed for E.Coli. Also performs gunzip and does not keep the original
        zip file. Uses wget.
        Args:
            ftp_list (lst) : array of ftp links already found, no ability to
                            call without specific ftp links
        Returns:
            None
    """

    # rename ncbi accession ids to common names
    output_names = ['HM27_FASTA.fna.gz', 'HM46_FASTA.fna.gz', \
                    'HM65_FASTA.fna.gz', 'HM69_FASTA.fna.gz', \
                    'hm27_feature.txt.gz', 'hm46_feature.txt.gz', \
                    'hm65_feature.txt.gz', 'hm69_feature.txt.gz']

    os.mkdir(FASTA_DIR_PATH)
    os.chdir(FASTA_DIR_PATH)
    # does not execute through shell
    for ftp_link, output_name in zip(ftp_list, output_names):
        wget_command = ['wget', '-O', output_name, ftp_link]
        gunzip_command = ['gunzip', output_name]

        subprocess.run(wget_command)
        subprocess.run(gunzip_command)

        LOGGER.info("Grabbed {}".format(output_name))
    os.chdir('../')


def parse_seqio_fasta(fasta_list, log_file):
    """ parses the refseq fasta files for each respective strain, counts
        the number of contigs per SeqIO record, and counts the total
        number of base pairs for each contig from each respective genome
        Args:
            fasta_record_list [lst] : array of SeqIO objects
            log_file [file] : output file needed to write to; opened at main()
        Returns:
            None
    """

    assembly_name_list = ['HM27', 'HM46', 'HM65', 'HM69']

    for fasta_file, assembly_name in zip(fasta_list, assembly_name_list):

        num_of_contigs = 0
        num_of_bp = 0

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            num_of_contigs += 1
            if len(seq_record.seq) >= 1000:
                num_of_bp += len(seq_record.seq)

        log_file.write('There are {} contigs in the {} assembly.\n'.format(\
                                                num_of_contigs, assembly_name))

        log_file.write('There are {} base pairs in the {} assembly.\n'.format(\
                                                    num_of_bp, assembly_name))


 #  _____  _____   ____  _  ___  __
 # |  __ \|  __ \ / __ \| |/ / |/ /    /\
 # | |__) | |__) | |  | | ' /| ' /    /  \
 # |  ___/|  _  /| |  | |  < |  <    / /\ \
 # | |    | | \ \| |__| | . \| . \  / ____ \
 # |_|    |_|  \_\\____/|_|\_\_|\_\/_/    \_\
 #

def build_prokka(fasta_list, log_file):
    """ takes a list of fasta paths of the four specified strains and
        runs the prokka software
        Utilizies the subprocess.run() method, but can be optimized to run
        subprocess.popen for parallel processing.
        Opens log file from main() and copies contents from the .txt file
        generated from the prokka output
        Args:
            fasta_list [lst] : array of fasta paths for prokka
        Returns:
            None
    """
    # prokka output directory names, do not change!
    prokka_output_dir = ['prokka_hm27', 'prokka_hm46', 'prokka_hm65', 'prokka_hm69']

    # base name for each prokka file in each prokka_output_directory
    genome_name_list = ['hm27_anno', 'hm46_anno', 'hm65_anno', 'hm69_anno']

    for fasta_file, output_dir, genome_name in zip(fasta_list, prokka_output_dir, genome_name_list):

        prokka_command = "prokka --outdir {} --prefix {} {} --genus Escherichia".format(output_dir, \
                                                                                        genome_name, \
                                                                                        fasta_file)

        LOGGER.info("Prokka Command Ran: {}".format(prokka_command))
        log_file.write("Prokka Command Ran: {}".format(prokka_command))
        subprocess.run(prokka_command, shell=True)

        LOGGER.info("Prokka_finished for {}".format(genome_name))

        txt_prokka_output = CURRENT_DIR + '/{}/{}.txt'.format(output_dir, genome_name)

        LOGGER.info("Output Directory Check {}".format(txt_prokka_output))
        LOGGER.info("Copying {}.txt contents to log file".format(genome_name))

        # copies prokka output text file to log file
        with open(txt_prokka_output, 'r') as txt_file:
            log_file.write("\n{} annotation\n".format(output_dir))
            for line in txt_file:
                log_file.write(str(line))




 #   _____ _____         _______ ____   ____  _       _____
 #  / ____|  __ \     /\|__   __/ __ \ / __ \| |     / ____|
 # | (___ | |__) |   /  \  | | | |  | | |  | | |    | (___
 #  \___ \|  _  /   / /\ \ | | | |  | | |  | | |     \___ \
 #  ____) | | \ \  / ____ \| | | |__| | |__| | |____ ____) |
 # |_____/|_|  \_\/_/    \_\_|  \____/ \____/|______|_____/
 #

def prefetch_fastq_decomp():
    """ utilizes tools from SRA toolkit to grab SRA files, and decompress them to
        fastq files. Need to improve with parallel decompression.
        Tools:
            prefetch : safely builds ftp path with sra file accession ids
                        grabs all required dependencies and reference files
                        and saves them to $HOME/ncbi/public/sra. DO NOT CHANGE.
            fastq-dump : decompresses the SRA files and saves them to directed
                        fastq directory 'id_sra'; searches the public directory
                        $HOME/ncbi/public/sra for these files, will redownload
                        if not present. DO NOT MOVE FILES.
        Args:
            None
        Returns:
            None
    """

    # SRA accession ids to feed to prefetch, will download all surrounding
    # dependency files.
    sra_files = ['SRR1278956', 'SRR1278960', 'SRR1283106', 'SRR1278963']

    # output directory names for fastq-dump; do not change
    fastq_dir_list = ['hm27_sra', 'hm46_sra', 'hm65_sra', 'hm69_sra']

    # prefetch places SRA files in $HOME/ncbi/public/sra
    # do not move
    for sraf in sra_files:
        LOGGER.info("Fetching {}".format(sraf))
        subprocess.run(['prefetch', sraf])
        LOGGER.info("Fetched {}".format(sraf))

    # fastq-dump looks in $HOME/ncbi/public/sra
    for sra_file, fastq_out in zip(sra_files, fastq_dir_list):

        LOGGER.info("Decompressing {}".format(sra_file))
        print("Decompressing {}".format(sra_file))

        fq_command = "fastq-dump -I --split-files {} -O {}".format(sra_file, fastq_out)
        subprocess.run(fq_command, shell=True)

        print("Decompression Complete")
        LOGGER.info("Finished {}".format(sra_file))



 #  _______ ____  _____  _    _       _______
 # |__   __/ __ \|  __ \| |  | |   /\|__   __|
 #    | | | |  | | |__) | |__| |  /  \  | |
 #    | | | |  | |  ___/|  __  | / /\ \ | |
 #    | | | |__| | |    | |  | |/ ____ \| |
 #    |_|  \____/|_|    |_|  |_/_/    \_\_|
 #

def build_tophat_alignment(fasta_file_list, gff_list, fastq_tuple_list, bam_file_list, sorted_bam_list, threads):
    """ main tophat/bowtie2 build. requires that the reference fasta files that were found
        are indexed by bowtie2. Once indexed by bowtie2, tophat will create a transcriptome
        index that is used for alignment. After that is built, tophat2 will perform the
        alignment ~ 3-4 hours. Once each alignment is performed, samtools is called to
        sort each bam in the same fashion and outputs file.sorted.bam file.
        Tools:
            bowtie2 : for indexing reference
            tophat2 : builds transcriptome index and performs RNA-seq alignment
            samtools : sorts the bams for downstream analysis
        Args:
            fasta_file_list (lst) : array of paths to fasta files
            gff_list (lst) : array of paths to gff files from prokka
            fastq_tuple_list (tup/lst) : holds the paths to fastq1 and fastq 2
                        from fastq dump
            bam_file_list (lst) : array of paths to where the bams will be located
                                    after alignment
            sorted_bam_list (lst) : array of path names for sorted bams after samtools
        Returns:
            None
    """

    idx_base_list = ['hm27_index', 'hm46_index', 'hm65_index', 'hm69_index']

    tophat_output_dir = ['hm27_tophat', 'hm46_tophat', 'hm65_tophat', 'hm69_tophat']

    trans_idx_list = ['transcriptome/hm27_index', 'transcriptome/hm46_index', \
    'transcriptome/hm65_index', 'transcriptome/hm69_index']
    # Begins to build the bowtie2 index for each reference sample
    # Must make a copy of the fasta file into the same format as base name
    # for tophat2, but with the .fa file type, NOT .fna.
    LOGGER.info("Beginning bowtie2 to build reference index.")
    for fna_file, base_name in zip(fasta_file_list, idx_base_list):

        bwt2_command = "bowtie2-build --threads {} -f {} {}".format(threads, fna_file, base_name)
        copy_command = "cp {} {}.fa".format(fna_file, base_name)
        subprocess.run(bwt2_command, shell=True)
        subprocess.run(copy_command, shell=True)
        LOGGER.info("Copied {} file to {}.fa".format(fna_file, base_name))
        LOGGER.info("Built reference index for {}".format(base_name))


    LOGGER.info("Beginning tophat to perform alignment.")
    # for gff_file, trans_idx, idx_base_name in zip(trans_idx_list, gff_list, idx_base_list):
    #
    #     trans_idx_command = "tophat -G {} --transcriptome-index={} {}".format(gff_file, \
    #                                                             trans_idx, idx_base_name)
    #     subprocess.run(trans_idx_command, shell=True)

    for tp_out_name, gff_file, idx_base_name, fastq_tup in zip(tophat_output_dir, gff_list, \
                                                                idx_base_list, fastq_tuple_list):

        top_hat_command = "tophat2 -p {} -o {} -G {} {} {} {}".format(threads, tp_out_name, \
                                                                        gff_file, idx_base_name, \
                                                                        fastq_tup[0], fastq_tup[1])

        LOGGER.info("Aligning {}".format(idx_base_name))

        subprocess.run(top_hat_command, shell=True)

        LOGGER.info("Alignment Complete")


    print("Sorting BAM files.")
    LOGGER.info("Beginning bam file sorting.")
    for bam_file, sorted_out_bam in zip(bam_file_list, sorted_bam_list):

        LOGGER.info("Sorting {}".format(bam_file))

        sort_bam_command = "samtools sort {} -o {}".format(bam_file, sorted_out_bam)
        subprocess.run(sort_bam_command, shell=True)

        LOGGER.info("Finished sorting {}".format(sorted_out_bam))



 #   _____ _    _ ______ ______ _      _____ _   _ _  __ _____
 #  / ____| |  | |  ____|  ____| |    |_   _| \ | | |/ // ____|
 # | |    | |  | | |__  | |__  | |      | | |  \| | ' /| (___
 # | |    | |  | |  __| |  __| | |      | | | . ` |  <  \___ \
 # | |____| |__| | |    | |    | |____ _| |_| |\  | . \ ____) |
 #  \_____|\____/|_|    |_|    |______|_____|_| \_|_|\_\_____/
 #

def run_cufflinks_suite(gff_list, sorted_bam_list, assembly_file, merged_gtf, threads):
    """ runs cufflinks suite to assemble the transcripts from the sorted
        bam files, merges them, and then normalizes the gene counts
        Tools:
            cufflinks
            cuffmerge
            cuffnorm
        Args:
            gff_list (lst) : array of paths to the gff files from prokka
            sorted_bam_list (lst) : array of paths to the sorted bams
            assembly_file
            merged_gtf (gtf) : merged from cuffmerge
        Returns:
            None
    """

    # array of top level directory for cufflink output files
    cuff_out_list = ['hm27_cuff', 'hm46_cuff', 'hm65_cuff', 'hm69_cuff']

    # building transcript assemblies with the sorted bam
    for gff_file, cuff_out, sorted_bam in zip(gff_list, cuff_out_list, sorted_bam_list):
        LOGGER.info("Assembling transcript for {}".format(cuff_out))

        cufflink_command = "cufflinks -p {} -G {} -o {} {}".format(threads, \
                                                                    gff_file, \
                                                                    cuff_out, \
                                                                    sorted_bam)

        subprocess.run(cufflink_command, shell=True)
        LOGGER.info("Assembly complete.")


    # merging each assembled transcriptome into 1 transcript
    LOGGER.info("Merging all assembled transcripts from each genome.")
    cuffmerge_command = "cuffmerge -p {} -o {} {}".format(threads, 'merged_ecoli', assembly_file)
    subprocess.run(cuffmerge_command, shell=True)
    LOGGER.info("Merging complete.")


    # normalizing the merged transcriptome against each of its sorted bam
    LOGGER.info("Normalizing the merged transcriptome.")
    cuffnorm_command = "cuffnorm -o diff_results -p {} {} {} {} {} {}".format(threads, merged_gtf, \
                                                                    sorted_bam_list[0], \
                                                                    sorted_bam_list[1], \
                                                                    sorted_bam_list[2], \
                                                                    sorted_bam_list[3])

    subprocess.run(cuffnorm_command, shell=True)
    LOGGER.info("Normalization complete.")



 #  __  __          _____ _   _
 # |  \/  |   /\   |_   _| \ | |
 # | \  / |  /  \    | | |  \| |
 # | |\/| | / /\ \   | | | . ` |
 # | |  | |/ ____ \ _| |_| |\  |
 # |_|  |_/_/    \_\_____|_| \_|
 #

def main():
    """ runs the main script in linear order, need to optimize for parallel processing"""

    args = arg_parser()
    threads = args.threads
    LOGGER.info("{} threads flagged to run each software.".format(threads))

    LOGGER.info("Pipeline beginning.")
    log_file = open('UPEC.log', 'w')

    # Storing all constants in lists for efficient looping within large functions
    # maintains the same order, and ensures no missing files

    # input for ncbi wget, holds ftp paths
    fasta_ftp_list = [HM27_FILES[0], HM46_FILES[0], HM65_FILES[0], HM69_FILES[0], \
                    HM27_FILES[1], HM46_FILES[1], HM65_FILES[1], HM69_FILES[1]]

    # output from ncbi wget
    fasta_file_list = [HM27_FASTA, HM46_FASTA, HM65_FASTA, HM69_FASTA]

    # output from fastq-dump
    fastq_tuple_list = [(HM27_FASTQ_1, HM27_FASTQ_2), \
                        (HM46_FASTQ_1, HM46_FASTQ_2), \
                        (HM65_FASTQ_1, HM65_FASTQ_2), \
                        (HM69_FASTQ_1, HM69_FASTQ_2)]

    # output from prokka software
    gff_list = [HM27_GFF_FILE, HM46_GFF_FILE, HM65_GFF_FILE, HM69_GFF_FILE]

    # output from tophat2
    bam_file_list = [HM27_BAM, HM46_BAM, HM65_BAM, HM69_BAM]

    # output from samtools sort and tophat2
    sorted_bam_list = [HM27_SORTED_BAM, HM46_SORTED_BAM, HM65_SORTED_BAM, HM69_SORTED_BAM]

    # grabbing ftp files from ncbi using wget method
    # LOGGER.info("Beginning to find FTP files.")
    # wget_gunzip_fasta(fasta_ftp_list)
    #
    # # Parsing the FASTA and counting number of contigs and base pairs > 1000 in length
    # LOGGER.info("Parsing FASTA and writing to log file.")
    # parse_seqio_fasta(fasta_file_list, log_file)
    #
    # LOGGER.info("Starting gene annotation with Prokka")
    # build_prokka(fasta_file_list, log_file)

    # LOGGER.info("Grabbing SRA files and converting to FASTQ")
    # prefetch_fastq_decomp()

    LOGGER.info("Beginning alignment process.")
    build_tophat_alignment(fasta_file_list, gff_list, fastq_tuple_list, \
                                                bam_file_list, sorted_bam_list, threads)

    LOGGER.info("Creating assembly file for cuffmerge.")
    with open('ecoli_assemblies.txt', 'w') as assemble_file:
        assemble_file.write("./hm27_cuff/transcripts.gtf\n")
        assemble_file.write("./hm46_cuff/transcripts.gtf\n")
        assemble_file.write("./hm65_cuff/transcripts.gtf\n")
        assemble_file.write("./hm69_cuff/transcripts.gtf\n")

    LOGGER.info("Beginning to run cufflinks")
    run_cufflinks_suite(gff_list, sorted_bam_list, 'ecoli_assemblies.txt', MERGED_GTF, threads)

    LOGGER.info("Pipeline Complete.")
    log_file.close()

if __name__ == '__main__':
    main()
