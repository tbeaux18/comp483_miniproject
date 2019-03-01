#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0


tophat2.py


"""


import argparse
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-a', '--assembly_name', help='base name of the strain file from prokka.')
    parser.add_argument('-f', '--fasta_file', help='path to fasta file.')
    parser.add_argument('-s', '--sra_file', help='sra accession id')
    parser.add_argument('-t', '--threads', help='number of threads to use')

    return parser.parse_args()



def build_tophat_alignment(fasta_file, idx_name, fastq1, fastq2, output_dir, threads):
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

    with open(sra_file, 'r') as sra_ids:
        for line in sra_ids:
            new_line = line.strip().split(',')
            assembly_name = new_line[0]
            sra_acc_id = new_line[1]
            fasta_file = './ncbi_fasta/' + assembly_name + '_FASTA.fna'
            output_dir = assembly_name + '_tophat'
            idx_name = assembly_name + '_index'
            fastq_1 = './' + assembly_name + '_sra/' + sra_acc_id + '_1.fastq'
            fastq_2 = './' + assembly_name + '_sra/' + sra_acc_id + '_2.fastq'

            build_tophat_alignment(fasta_file, idx_name, fastq_1, fastq_2, output_dir, threads)

def main():
    args = arg_parser()

    sra_file = args.sra_file
    threads = args.threads

    build_tophat_run(sra_file, threads)

if __name__ == '__main__':
    main()
