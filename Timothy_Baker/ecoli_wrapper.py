#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Timothy Baker
@version: 1.0.0

ecoli_tb.py

This module is a python wrapper that pulls the FASTA files from specific FTP URLs from NCBI.

ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz HM27
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz HM46
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz HM65
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz HM69

ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz HM27
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz HM46
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz HM65
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz HM69

2. Calculate the number of contigs for each assembly and write the # out to the log file as follows:
There are # contigs in the assembly HM27.
[You’ll likewise write out the number of contigs for the other 3 strains.]

3. Calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length) and write this # out to the log file as follows:
There are # bp in the assembly HM27.
[You’ll likewise write out the number of bp for the other 3 strains.]

4. Use Prokka to annotate these assembly. Since we’re working on an E. coli genome, we are in luck… there’s already an Escherichia genus database. Let’s use it. Write the Prokka command to the log file.

5. Write the results of the annotation in the *.txt file to the log file in the same format as the *.txt file for each strain; add a label so someone reading the log could see which results correspond to which strain.

6. The assembled genome in RefSeq for E. coli K-12 (NC_000913) has 4140 CDS and 89 tRNAs annotated. Write to the log file the discrepancy (if any) found. For instance, if my Prokka annotation predicted 4315 CDS and 88 tRNA’s, I would write,
Prokka found 175 additional CDS and 1 less tRNA than the RefSeq in assembly HM27.
[You’ll likewise write out the number of bp for the other 3 strains.]

7. Now that we know where the genes are located, we can see how these genes are transcribed.
Use TopHat & Cufflinks to map the reads of a specific strain to the genome of the strain and quantify their expression, respectively.
Details re: the processes of TopHat and Cufflinks can be found in the class slides and Trapnell et al. 2013 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/.
 The accession numbers for the transcriptomes are listed below,
HM27:  https://www.ncbi.nlm.nih.gov/sra/SRX541301
ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278956/SRR1278956.sra
HM46:  https://www.ncbi.nlm.nih.gov/sra/SRX541306
ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra
HM65:  https://www.ncbi.nlm.nih.gov/sra/SRX541312
ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra
HM69:  https://www.ncbi.nlm.nih.gov/sra/SRX541316
ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra

9. Using Cuffnorm, normalize all 4 of your transcriptomes. Parse the output of Cuffnorm (the Simple-table gene attributes format file) such that you create a sorted file, with the highest expressed gene first, for each strain and write this to file, e.g. HM27_normalized_sorted.tsv.

Dependencies:
    Prokka
    tophat2
    bowtie2
    sratools
        prefetch
        fastq-dump
    
    Biopython

"""

import os
import subprocess
import logging
from Bio import SeqIO

CURRENT_DIR = os.getcwd()

HM27_FASTA = 'HM27_FASTA.fna'
HM46_FASTA = 'HM46_FASTA.fna'
HM65_FASTA = 'HM65_FASTA.fna'
HM69_FASTA = 'HM69_FASTA.fna'

HM27_GFF_FILE = CURRENT_DIR + '/prokka_hm27/hm27_index.gff'
HM46_GFF_FILE = CURRENT_DIR + '/prokka_hm46/hm46_index.gff'
HM65_GFF_FILE = CURRENT_DIR + '/prokka_hm65/hm65_index.gff'
HM69_GFF_FILE = CURRENT_DIR + '/prokka_hm69/hm69_index.gff'

HM27_BAM = CURRENT_DIR + '/hm27_tophat/accepted_hits.bam'
HM46_BAM = CURRENT_DIR + '/hm46_tophat/accepted_hits.bam'
HM65_BAM = CURRENT_DIR + '/hm65_tophat/accepted_hits.bam'
HM69_BAM = CURRENT_DIR + '/hm69_tophat/accepted_hits.bam'

HM27_SORTED_BAM = CURRENT_DIR + '/hm27_tophat/accepted_hits.sorted.bam'
HM46_SORTED_BAM = CURRENT_DIR + '/hm46_tophat/accepted_hits.sorted.bam'
HM65_SORTED_BAM = CURRENT_DIR + '/hm65_tophat/accepted_hits.sorted.bam'
HM69_SORTED_BAM = CURRENT_DIR + '/hm69_tophat/accepted_hits.sorted.bam'

HM27_FASTQ_1 = CURRENT_DIR + '/hm27_sra/SRR1278956_1.fastq'
HM27_FASTQ_2 = CURRENT_DIR + '/hm27_sra/SRR1278956_2.fastq'
HM46_FASTQ_1 = CURRENT_DIR + '/hm46_sra/SRR1278960_1.fastq'
HM46_FASTQ_2 = CURRENT_DIR + '/hm46_sra/SRR1278960_2.fastq'
HM65_FASTQ_1 = CURRENT_DIR + '/hm65_sra/SRR1283106_1.fastq'
HM65_FASTQ_2 = CURRENT_DIR + '/hm65_sra/SRR1283106_2.fastq'
HM69_FASTQ_1 = CURRENT_DIR + '/hm69_sra/SRR1278963_1.fastq'
HM69_FASTQ_2 = CURRENT_DIR + '/hm69_sra/SRR1278963_2.fastq'

# FTP FILES NEEDED FOR ANALYSIS, TUPLE FORMAT, (FASTA FTP LINK, FEATURE COUNT FTP LINK)
HM27_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz')

HM46_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz')

HM65_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz')

HM69_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz')



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

    # fixed output names
    output_names = ['HM27_FASTA.fna.gz', 'HM46_FASTA.fna.gz', \
                    'HM65_FASTA.fna.gz', 'HM69_FASTA.fna.gz', \
                    'hm27_feature.txt.gz', 'hm46_feature.txt.gz', \
                    'hm65_feature.txt.gz', 'hm69_feature.txt.gz']

    # does not execute through shell
    for ftp_link, output_name in zip(ftp_list, output_names):
        wget_command = ['wget', '-O', output_name, ftp_link]
        gunzip_command = ['gunzip', output_name]

        subprocess.run(wget_command)

        subprocess.run(gunzip_command)

        LOGGER.info("Grabbed {}".format(output_name))



def parse_seqio_fasta(fasta_record_list, log_file):
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

    for single_fasta_record, assembly_name in zip(fasta_record_list, \
                                                        assembly_name_list):
        num_of_contigs = len(single_fasta_record)

        num_of_bp = 0

        for individ_record in single_fasta_record:

            num_of_bp += len(individ_record.seq)

        log_file.write('There are {} contigs in the {} assembly.\n'.format(\
                                                num_of_contigs, assembly_name))

        log_file.write('There are {} bp in the {} assembly.\n'.format(\
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
        log_file.write("Prokka Command Ran: s{}".format(prokka_command))
        subprocess.run(prokka_command, shell=True)

        LOGGER.info("Prokka_finished for {}".format(genome_name))

        txt_prokka_output = CURRENT_DIR + '/{}/{}.txt'.format(output_dir, genome_name)

        LOGGER.info("Output Directory Check {}".format(txt_prokka_output))
        LOGGER.info("Copying {}.txt contents to log file".format(genome_name))

        # copies prokka output text file to log file
        with open(txt_prokka_output, 'r') as txt_file:
            with open(log_file, "w") as output_file:
                output_file.write("\n{} annotation\n".format(output_dir))
                for line in txt_file:
                    output_file.write(line)




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

    fastq_dir_list = ['hm27_sra', 'hm46_sra', 'hm65_sra', 'hm69_sra']

    # prefetch places SRA files in $HOME/ncbi/public/sra
    # do not move
    for sraf in sra_files:
        LOGGER.info("Fetching {}".format(sraf))
        subprocess.run(['prefetch', sraf])
        LOGGER.info("Fetched {}".format(sraf))

    for sra_file, fastq_out in zip(sra_files, fastq_dir_list):

        LOGGER.info("Decompressing {}".format(sra_file))

        fq_command = "fastq-dump -I --split-files {} -O {}".format(sra_file, fastq_out)

        print("Decompressing {}".format(sra_file))

        # fastq-dump looks in $HOME/ncbi/public/sra
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

def build_tophat_alignment(fasta_file_list, gff_list, fastq_tuple_list, bam_file_list, sorted_bam_list):
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

    transcriptome_idx = ['transcriptome/hm27', 'transcriptome/hm46', \
                        'transcriptome/hm65', 'transcriptome/hm69']

    # Begins to build the bowtie2 index for each reference sample
    # Must make a copy of the fasta file into the same format as base name
    # for tophat2, but with the .fa file type, NOT .fna.
    LOGGER.info("Beginning bowtie2 to build reference index.")
    for fna_file, base_name in zip(fasta_file_list, idx_base_list):

        bwt2_command = "bowtie2-build --threads 2 -f {} {}".format(fna_file, base_name)
        copy_command = "cp fna_file {}.fa".format(base_name)
        subprocess.run(bwt2_command, shell=True)
        subprocess.run(copy_command, shell=True)
        LOGGER.info("Copied {} file to {}.fa".format(fna_file, base_name))
        LOGGER.info("Built reference index for {}".format(base_name))

    LOGGER.info("Beginning tophat to perform alignment.")
    for gff_file, idx_base_name, trans_idx, tp_out_name, fastq_tup in zip(gff_list, \
                                                            idx_base_list, \
                                                            transcriptome_idx, \
                                                            tophat_output_dir, \
                                                            fastq_tuple_list):

        trans_idx_command = "tophat -G {} --transcriptome-index={} {}".format(gff_file, \
                                                            trans_idx, \
                                                            idx_base_name)

        top_hat_command = "tophat2 -p 4 --transcriptome-index={} {} -o {} {} {} {}".format(trans_idx, idx_base_name, \
                                                                                            tp_out_name, idx_base_name, \
                                                                                        fastq_tup[0], fastq_tup[1])

        LOGGER.info("Aligning {}".format(idx_base_name))

        subprocess.run(trans_idx_command, shell=True)
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

def run_cufflinks_suite(gff_list, sorted_bam_list, assembly_file, merged_gtf):
    cuff_out_list = ['hm27_cuff', 'hm46_cuff', 'hm65_cuff', 'hm69_cuff']

    for gff_file, cuff_out, sorted_bam in zip(gff_list, cuff_out_list, sorted_bam_list):
        cufflink_command = "cufflinks -p 4 -G {} -o {} {}".format(gff_file, cuff_out, sorted_bam)
        subprocess.run(cufflink_command, shell=True)

    cuffmerge_command = "cuffmerge -p 4 -o {} {}".format('merged_ecoli', assembly_file)
    subprocess.run(cuffmerge_command, shell=True)

    cuffnorm_command = "cuffnorm -o diff_results -p 4 {} {} {} {} {}".format(merged_gtf, \
                                                                    sorted_bam_list[0], \
                                                                    sorted_bam_list[1], \
                                                                    sorted_bam_list[2], \
                                                                    sorted_bam_list[3])

    subprocess.run(cuffnorm_command, shell=True)



 #  __  __          _____ _   _
 # |  \/  |   /\   |_   _| \ | |
 # | \  / |  /  \    | | |  \| |
 # | |\/| | / /\ \   | | | . ` |
 # | |  | |/ ____ \ _| |_| |\  |
 # |_|  |_/_/    \_\_____|_| \_|
 #

def main():

    log_file = open('UPEC.log', 'w')
    # when script runs, need to set the current working directory
    # create the environment variables first of the what the system needs to do
    # create the appropriate directory on the new system
    # grab first the working directory so you can properly manage where
    # each of the files go
    # cwd = os.getcwd()

    # hm27_fasta, hm46_fasta, hm65_fasta, hm69_fasta = 'HM27_FASTA.fna', \
    #                                                     'HM46_FASTA.fna', \
    #                                                     'HM65_FASTA.fna', \
    #                                                     'HM69_FASTA.fna'
    #
    # hm27_gff_file = cwd + '/prokka_hm27/hm27_index.gff'
    # hm46_gff_file = cwd + '/prokka_hm46/hm46_index.gff'
    # hm65_gff_file = cwd + '/prokka_hm65/hm65_index.gff'
    # hm69_gff_file = cwd + '/prokka_hm69/hm69_index.gff'
    # hm27_bam = cwd + '/hm27_tophat/accepted_hits.bam'
    # hm46_bam = cwd + '/hm46_tophat/accepted_hits.bam'
    # hm65_bam = cwd + '/hm65_tophat/accepted_hits.bam'
    # hm69_bam = cwd + '/hm69_tophat/accepted_hits.bam'
    # hm27_sorted_bam = cwd + '/hm27_tophat/accepted_hits.sorted.bam'
    # hm46_sorted_bam = cwd + '/hm46_tophat/accepted_hits.sorted.bam'
    # hm65_sorted_bam = cwd + '/hm65_tophat/accepted_hits.sorted.bam'
    # hm69_sorted_bam = cwd + '/hm69_tophat/accepted_hits.sorted.bam'
    # hm27_fastq_1 = cwd + '/hm27_sra/SRR1278956_1.fastq'
    # hm27_fastq_2 = cwd + '/hm27_sra/SRR1278956_2.fastq'
    # hm46_fastq_1 = cwd + '/hm46_sra/SRR1278960_1.fastq'
    # hm46_fastq_2 = cwd + '/hm46_sra/SRR1278960_2.fastq'
    # hm65_fastq_1 = cwd + '/hm65_sra/SRR1283106_1.fastq'
    # hm65_fastq_2 = cwd + '/hm65_sra/SRR1283106_2.fastq'
    # hm69_fastq_1 = cwd + '/hm69_sra/SRR1278963_1.fastq'
    # hm69_fastq_2 = cwd + '/hm69_sra/SRR1278963_2.fastq'

    fastq_tuple_list = [(HM27_FASTQ_1, HM27_FASTQ_2), \
                        (HM46_FASTQ_1, HM46_FASTQ_2), \
                        (HM65_FASTQ_1, HM65_FASTQ_2), \
                        (HM69_FASTQ_1, HM69_FASTQ_2)]

    gff_list = [HM27_GFF_FILE, HM46_GFF_FILE, HM65_GFF_FILE, HM69_GFF_FILE]
    fasta_file_list = [HM27_FASTA, HM46_FASTA, HM65_FASTA, HM69_FASTA]
    bam_file_list = [HM27_BAM, HM46_BAM, HM65_BAM, HM69_BAM]
    sorted_bam_list = [HM27_SORTED_BAM, HM46_SORTED_BAM, HM65_SORTED_BAM, HM69_SORTED_BAM]

    merged_gtf = CURRENT_DIR + '/merged_ecoli/merged.gtf'

    fasta_ftp_list = [HM27_FILES[0], HM46_FILES[0], HM65_FILES[0], HM69_FILES[0], \
                    HM27_FILES[1], HM46_FILES[1], HM65_FILES[1], HM69_FILES[1]]


    # # grabbing ftp files
    # LOGGER.info("Beginning to find FTP files.")
    # wget_gunzip_fasta(fasta_ftp_list)
    #
    # LOGGER.info("Starting gene annotation with Prokka")
    # build_prokka(fasta_file_list)
    #
    # LOGGER.info("Grabbing SRA files and converting to FASTQ")
    # prefetch_fastq_decomp()

    # need to add copy commmand to make the fasta files the same base name as bwt base
    # create alternative directory structure to consider this
    # need to configure to grab fastq files from specific directories and
    # store those paths as simple variables to pass through.
    # may need to change .gff name to the same base name, much of top hats functionality
    # is not kept up

    # LOGGER.info("Beginning alignment process.")
    # build_tophat_alignment(fasta_file_list, gff_list, fastq_tuple_list, \
    #                                             bam_file_list, sorted_bam_list)

    # with open('ecoli_assemblies.txt', 'w') as assemble:
    #     assemble.write("./hm27_cuff/transcripts.gtf\n./hm46_cuff/transcripts.gtf\n./hm65_cuff/transcripts.gtf\n/hm69_cuff/transcripts.gtf\n")
    #
    # LOGGER.info("Beginning to run cufflinks")
    # run_cufflinks_suite(gff_list, sorted_bam_list, 'ecoli_assemblies.txt', merged_gtf)

    # need to include grabbing file path names
    # think of how to store these records
    hm27_records = list(SeqIO.parse("HM27_FASTA.fna", "fasta"))
    hm46_records = list(SeqIO.parse("HM46_FASTA.fna", "fasta"))
    hm65_records = list(SeqIO.parse("HM65_FASTA.fna", "fasta"))
    hm69_records = list(SeqIO.parse("HM69_FASTA.fna", "fasta"))


    # Store all FASTA records in a list for quick retrieval and looping
    # fasta_record_list = [hm27_records, hm46_records, hm65_records, hm69_records]
    # parse_seqio_fasta(fasta_record_list, fasta_record_output, log_file)

    # log_file.write("\nHM27 Prokka Annotation\n")
    # subprocess.run("cat hm27-prokka-output.txt >> UPEC.log", shell=True)


    log_file.close()

if __name__ == '__main__':
    main()
