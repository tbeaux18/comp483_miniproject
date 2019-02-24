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
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz')

HM46_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz')

HM65_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz')

HM69_FILES = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz', \
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz')



LOGGER = logging.getLogger(__name__)

LOGGER.setLevel(logging.INFO)

FORMATTER = logging.Formatter('%(levelname)s:%(name)s:%(asctime)s:%(message)s')

FILE_HANDLER = logging.FileHandler("working.log")

FILE_HANDLER.setFormatter(FORMATTER)

LOGGER.addHandler(FILE_HANDLER)

# ask for number of threads


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



def build_prokka(fasta_list):

    prokka_output_dir = ['prokka_hm27', 'prokka_hm46', 'prokka_hm65', 'prokka_hm69']

    genome_name_list = ['hm27_anno', 'hm46_anno', 'hm65_anno', 'hm69_anno']

    for fasta_file, output_dir, genome_name in zip(fasta_list, prokka_output_dir, genome_name_list):

        prokka_command = "prokka --outdir {} --prefix {} {} --genus Escherichia".format(output_dir, genome_name, fasta_file)

        LOGGER.info("Prokka Command Ran: {}".format(prokka_command))

        subprocess.run(prokka_command, shell=True)




def fastq_decomp():

    # local variables, do not overwrite, hardcoded directories for further
    # downstream analysis
    sra_files = ['SRR1278956', 'SRR1278960', 'SRR1283106', 'SRR1278963']

    fastq_dir_list = ['hm27_sra', 'hm46_sra', 'hm65_sra', 'hm69_sra']

    # prefetch places SRA files in $HOME/ncbi/public/sra
    # do not move
    for sraf in sra_files:
        subprocess.run(['prefetch', sraf])


    for sra_file, fastq_out in zip(sra_files, fastq_dir_list):

        LOGGER.info("Decompressing {}".format(sra_file))

        fq_command = "fastq-dump -I --split-files {} -O {}".format(sra_file, fastq_out)

        print("Decompressing {}".format(sra_file))
        subprocess.run(fq_command, shell=True)

        print("Decompression Complete")
        LOGGER.info("Finished {}".format(sra_file))





def wget_gunzip_fasta(ftp_list, output_list):

    fasta_output_name = ['HM27_FASTA.fna.gz', 'HM46_FASTA.fna.gz', \
                            'HM65_FASTA.fna.gz', 'HM69_FASTA.fna.gz']

    feature_txt_output = ['hm27_feature.txt.gz', 'hm46_feature.txt.gz', \
                            'hm65_feature.txt.gz', 'hm69_feature.txt.gz']

    for ftp_link, output_name in zip(ftp_list, output_list):
        wget_command = ['wget', '-O', output_name, ftp_link]

        gunzip_command = ['gunzip', output_name]

        subprocess.run(wget_command)

        subprocess.run(gunzip_command)





def build_tophat_alignment(fasta_file_list, gff_list, fastq_tuple_list, bam_file_list, sorted_bam_list):

    idx_base_list = ['hm27_index', 'hm46_index', 'hm65_index', 'hm69_index']

    tophat_output_dir = ['hm27_tophat', 'hm46_tophat', 'hm65_tophat', 'hm69_tophat']

    # Begins to build the bowtie2 index for each reference sample
    for fna_file, base_name in zip(fasta_file_list, idx_base_list):

        bwt2_command = "bowtie2-build --threads 2 -f {} {}".format(fna_file, base_name)

        subprocess.run(bwt2_command, shell=True)


    for gff_file, idx_base_name, tp_out_name, fastq_tup in zip(gff_list, \
                                                            idx_base_list, \
                                                            tophat_output_dir, \
                                                            fastq_tuple_list):

        trans_idx_command = "tophat -G {} --transcriptome-index={} {}".format(gff_file, \
                                                            idx_base_name, \
                                                            idx_base_name)

        top_hat_command = "tophat2 -p 4 -o {} {} {} {}".format(tp_out_name, \
                                                            idx_base_name, \
                                                            fastq_tup[0], \
                                                            fastq_tup[1])

        LOGGER.info("Aligning {}".format(idx_base_name))

        subprocess.run(trans_idx_command, shell=True)
        subprocess.run(top_hat_command, shell=True)

        LOGGER.info("Alignment Complete")

    for bam_file, sorted_out_bam in zip(bam_file_list, sorted_bam_list):

        sort_bam_command = "samtools sort {} -o {}".format(bam_file, sorted_out_bam)
        subprocess.run(sort_bam_command, shell=True)




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



def main():

    log_file = open('UPEC.log', 'w')
    # when script runs, need to set the current working directory
    # create the environment variables first of the what the system needs to do
    # create the appropriate directory on the new system
    # grab first the working directory so you can properly manage where
    # each of the files go
    cwd = os.getcwd()

    hm27_fasta, hm46_fasta, hm65_fasta, hm69_fasta = 'HM27_FASTA.fna', \
                                                        'HM46_FASTA.fna', \
                                                        'HM65_FASTA.fna', \
                                                        'HM69_FASTA.fna'

    hm27_filename = 'HM27_FASTA.fna'
    hm46_filename = 'HM46_FASTA.fna'
    hm65_filename = 'HM65_FASTA.fna'
    hm69_filename = 'HM69_FASTA.fna'
    hm27_gff_file = cwd + '/prokka_hm27/hm27_index.gff'
    hm46_gff_file = cwd + '/prokka_hm46/hm46_index.gff'
    hm65_gff_file = cwd + '/prokka_hm65/hm65_index.gff'
    hm69_gff_file = cwd + '/prokka_hm69/hm69_index.gff'
    hm27_bam = cwd + '/hm27_tophat/accepted_hits.bam'
    hm46_bam = cwd + '/hm46_tophat/accepted_hits.bam'
    hm65_bam = cwd + '/hm65_tophat/accepted_hits.bam'
    hm69_bam = cwd + '/hm69_tophat/accepted_hits.bam'
    hm27_sorted_bam = cwd + '/hm27_tophat/accepted_hits.sorted.bam'
    hm46_sorted_bam = cwd + '/hm46_tophat/accepted_hits.sorted.bam'
    hm65_sorted_bam = cwd + '/hm65_tophat/accepted_hits.sorted.bam'
    hm69_sorted_bam = cwd + '/hm69_tophat/accepted_hits.sorted.bam'
    hm27_fastq_1 = cwd + '/hm27_sra/SRR1278956_1.fastq'
    hm27_fastq_2 = cwd + '/hm27_sra/SRR1278956_2.fastq'
    hm46_fastq_1 = cwd + '/hm46_sra/SRR1278960_1.fastq'
    hm46_fastq_2 = cwd + '/hm46_sra/SRR1278960_2.fastq'
    hm65_fastq_1 = cwd + '/hm65_sra/SRR1283106_1.fastq'
    hm65_fastq_2 = cwd + '/hm65_sra/SRR1283106_2.fastq'
    hm69_fastq_1 = cwd + '/hm69_sra/SRR1278963_1.fastq'
    hm69_fastq_2 = cwd + '/hm69_sra/SRR1278963_2.fastq'

    fastq_tuple_list = [(hm27_fastq_1, hm27_fastq_2), \
                        (hm46_fastq_1, hm46_fastq_2), \
                        (hm65_fastq_1, hm65_fastq_2), \
                        (hm69_fastq_1, hm69_fastq_2)]

    gff_list = [hm27_gff_file, hm46_gff_file, hm65_gff_file, hm69_gff_file]
    fasta_file_list = [hm27_fasta, hm46_fasta, hm65_fasta, hm69_fasta]
    bam_file_list = [hm27_bam, hm46_bam, hm65_bam, hm69_bam]
    sorted_bam_list = [hm27_sorted_bam, hm46_sorted_bam, hm65_sorted_bam, hm69_sorted_bam]

    merged_gtf = cwd + '/merged_ecoli/merged.gtf'

    # fasta_ftp_list = [HM27_FILES[0], HM46_FILES[0], HM65_FILES[0], HM69_FILES[0]]
    # wget_gunzip_fasta(fasta_ftp_list, fasta_output_name)


    # feature_ftp_list = [HM27_FILES[1], HM46_FILES[1], HM65_FILES[1], HM69_FILES[1]]

    # wget_gunzip_fasta(feature_ftp_list, feature_txt_output)
    #


    # need to add copy commmand to make the fasta files the same base name as bwt base
    # create alternative directory structure to consider this
    # need to configure to grab fastq files from specific directories and
    # store those paths as simple variables to pass through.
    # may need to change .gff name to the same base name, much of top hats functionality
    # is not kept up



    with open('ecoli_assemblies.txt', 'w') as assemble:
        assemble.write("./hm27_cuff/transcripts.gtf\n./hm46_cuff/transcripts.gtf\n./hm65_cuff/transcripts.gtf\n/hm69_cuff/transcripts.gtf\n")


    # need to include grabbing file path names
    # think of how to store these records
    # hm27_records = list(SeqIO.parse("HM27_FASTA.fna", "fasta"))
    # hm46_records = list(SeqIO.parse("HM46_FASTA.fna", "fasta"))
    # hm65_records = list(SeqIO.parse("HM65_FASTA.fna", "fasta"))
    # hm69_records = list(SeqIO.parse("HM69_FASTA.fna", "fasta"))
    #
    # #
    # # Store all FASTA records in a list for quick retrieval and looping
    # fasta_record_list = [hm27_records, hm46_records, hm65_records, hm69_records]
    # fasta_record_output = ['HM27', 'HM46', 'HM65', 'HM69']
    # parse_seqio_fasta(fasta_record_list, fasta_record_output, log_file)
    #
    # log_file.write("\nHM27 Prokka Annotation\n")
    # subprocess.run("cat hm27-prokka-output.txt >> UPEC.log", shell=True)


    log_file.close()

if __name__ == '__main__':
    main()
