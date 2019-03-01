#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

tophat2.py

"""

 #  _______ ____  _____  _    _       _______
 # |__   __/ __ \|  __ \| |  | |   /\|__   __|
 #    | | | |  | | |__) | |__| |  /  \  | |
 #    | | | |  | |  ___/|  __  | / /\ \ | |
 #    | | | |__| | |    | |  | |/ ____ \| |
 #    |_|  \____/|_|    |_|  |_/_/    \_\_|
 #


import argparse
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-s', '--sra_file', help='sra accession id')
    parser.add_argument('-t', '--threads', help='number of threads to use')
    return parser.parse_args()



def build_tophat_alignment(fasta_file, idx_name, fastq1, fastq2, output_dir, threads):
    """ builds tophat alignment and runs the command. Only works for PE reads
        Tools:
            bowtie2 : for indexing reference
            tophat2 : builds transcriptome index and performs RNA-seq alignment
            samtools : sorts the bams for downstream analysis
        Args:
            fasta_file (str) : path to fasta
            idx_name (str) : base name
            fastq1 (str) : path to fastq1
            fastq2 (str) : path to fastq2
            output_dir (str) : output directory
            threads (int) : number of threads to run on
        Returns:
            None
    """

    # Begins to build the bowtie2 index for each reference sample
    # Must make a copy of the fasta file into the same format as base name
    # for tophat2, but with the .fa file type, NOT .fna.
    bwt2_command = "bowtie2-build --threads {} -f {} {}".format(threads, \
                                                                    fasta_file, \
                                                                    idx_name)

    copy_fasta = "cp {} {}.fa".format(fasta_file, idx_name)

    subprocess.run(bwt2_command, shell=True)

    subprocess.run(copy_fasta, shell=True)

    top_hat_command = "tophat2 -p {} -o {} {} {} {}".format(threads, \
                                                                output_dir, \
                                                                idx_name, \
                                                                fastq1, \
                                                                fastq2)


    subprocess.run(top_hat_command, shell=True)


def build_tophat_run(sra_file, threads):
    """ builds the top hat run and builds the paths
        Args:
            sra_file (str) : input file
            threads (int) : threads to use
        Returns:
            None
    """

    with open(sra_file, 'r') as sra_ids:
        for line in sra_ids:

            # creates a list
            new_line = line.strip().split(',')

            # strain name is first; specified in input instruction
            assembly_name = new_line[0]

            # requires the SRA file accession ID; usually begins with SSR
            sra_acc_id = new_line[1]

            # builds the fasta_file path name
            fasta_file = './ncbi_fasta/' + assembly_name + '_FASTA.fna'

            # builds the output directory name using strain name
            output_dir = assembly_name + '_tophat'

            # builds the index base name for bowtie2
            idx_name = assembly_name + '_index'

            # builds the fastq files that are located in the _sra direcotry
            fastq_1 = './' + assembly_name + '_sra/' + sra_acc_id + '_1.fastq'
            fastq_2 = './' + assembly_name + '_sra/' + sra_acc_id + '_2.fastq'

            # runs the software
            build_tophat_alignment(fasta_file, idx_name, fastq_1, fastq_2, output_dir, threads)


def main():
    """ runs the main script """

    # set args
    args = arg_parser()

    # sra file is needed in specified format
    sra_file = args.sra_file
    threads = args.threads

    # runs the build; no error handling available.
    # takes roughly 3-4 hours per genome
    build_tophat_run(sra_file, threads)

if __name__ == '__main__':
    main()
