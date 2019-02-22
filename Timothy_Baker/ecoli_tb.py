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

8. Using Cuffdiff identify significant changes in transcript expression between the four transcriptomes. We expect, since all 4 are UPEC strains that their expression will be similar, but maybe not. From the Cuffdiff output file, read the cds.diff file. Write to your log each row in which column “significant” is “yes”.

9. Using Cuffnorm, normalize all 4 of your transcriptomes. Parse the output of Cuffnorm (the Simple-table gene attributes format file) such that you create a sorted file, with the highest expressed gene first, for each strain and write this to file, e.g. HM27_normalized_sorted.tsv.

Dependencies:
    Prokka
    Biopython

"""

import os
import subprocess
import logging
from Bio import SeqIO





# FTP FILES NEEDED FOR ANALYSIS, TUPLE FORMAT, (FASTA FTP LINK, FEATURE COUNT FTP LINK, SRA FTP LINK)
HM27_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra')

HM46_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra')

HM65_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra')

HM69_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra')


LOGGER = logging.getLogger(__name__)

LOGGER.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')

FILE_HANDLER = logging.FileHandler("working.log")

FILE_HANDLER.setFormatter(FORMATTER)

LOGGER.addHandler(FILE_HANDLER)



def parse_seqio_fasta(fasta_record_list, assembly_name_list, log_file):

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



def build_prokka(fasta_list, output_dir_list, genome_name_list):

    for fasta_file, output_dir, genome_name in zip(fasta_list, output_dir_list, genome_name_list):

        prokka_command = ["prokka --outdir {} \
                            --prefix {}, {}, '--genus', \
                            'Escherichia".format(output_dir, genome_name, fasta_file)]

        subprocess.run(prokka_command, shell=True)



def wget_gunzip_fasta(ftp_list, output_list):

    for ftp_link, output_name in zip(ftp_list, output_list):
        wget_command = ['wget', '-O', output_name, ftp_link]

        gunzip_command = ['gunzip', output_name]

        subprocess.run(wget_command)

        subprocess.run(gunzip_command)



def main():

    log_file = open('UPEC.log', 'w')
    # when script runs, need to set the current working directory
    # create the environment variables first of the what the system needs to do
    # create the appropriate directory on the new system
    # grab first the working directory so you can properly manage where
    # each of the files go
    cwd = os.getcwd()



    hm27_filename = cwd + '/HM27_fasta.fna'
    hm46_filename = cwd + '/HM46_fasta.fna'
    hm65_filename = cwd + '/HM65_fasta.fna'
    hm69_filename = cwd + '/HM69_fasta.fna'


    # fasta_ftp_list = [HM27_FILES[0], HM46_FILES[0], HM65_FILES[0], HM69_FILES[0]]
    # fasta_output_name = ['HM27_FASTA.fna.gz', 'HM46_FASTA.fna.gz', \
    #                         'HM65_FASTA.fna.gz', 'HM69_FASTA.fna.gz']
    # wget_gunzip_fasta(fasta_ftp_list, fasta_output_name)

    filenames_list = [hm27_filename, hm46_filename, hm65_filename, hm69_filename]
    output_dir_list = ['prokka_hm27', 'prokka_hm46', 'prokka_hm65', 'prokka_hm69']
    genome_names = ['hm27_anno', 'hm46_anno', 'hm65_anno', 'hm69_anno']
    build_prokka(filenames_list, output_dir_list, genome_names)


    feature_ftp_list = [HM27_FILES[1], HM46_FILES[1], HM65_FILES[1], HM69_FILES[1]]
    feature_txt_output = ['hm27_feature.txt.gz', 'hm46_feature.txt.gz', \
                            'hm65_feature.txt.gz', 'hm69_feature.txt.gz']
    wget_gunzip_fasta(feature_ftp_list, feature_txt_output)




    # need to include grabbing file path names
    # think of how to store these records
    hm27_records = list(SeqIO.parse("HM27_FASTA.fna", "fasta"))
    hm46_records = list(SeqIO.parse("HM46_FASTA.fna", "fasta"))
    hm65_records = list(SeqIO.parse("HM65_FASTA.fna", "fasta"))
    hm69_records = list(SeqIO.parse("HM69_FASTA.fna", "fasta"))


    # Store all FASTA records in a list for quick retrieval and looping
    fasta_record_list = [hm27_records, hm46_records, hm65_records, hm69_records]
    fasta_record_output = ['HM27', 'HM46', 'HM65', 'HM69']
    parse_seqio_fasta(fasta_record_list, fasta_record_output, log_file)




    log_file.close()

if __name__ == '__main__':
    main()
