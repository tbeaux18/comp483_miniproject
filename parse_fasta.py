#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

"""

import os
import argparse
import subprocess
from Bio import SeqIO





def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-f', '--ftp_file', help='path to text file that has ftp links')

    return parser.parse_args()


 #  ______       _____ _______
 # |  ____/\    / ____|__   __|/\
 # | |__ /  \  | (___    | |  /  \
 # |  __/ /\ \  \___ \   | | / /\ \
 # | | / ____ \ ____) |  | |/ ____ \
 # |_|/_/    \_\_____/   |_/_/    \_\
 #



def wget_gunzip_fasta(output_name, ftp_link):
    """ grabs the zipped fasta and feature txt files for the specific strains
        listed for E.Coli. Also performs gunzip and does not keep the original
        zip file. Uses wget.
        Args:
            ftp_list (lst) : array of ftp links already found, no ability to
                            call without specific ftp links
        Returns:
            None
    """

    if os.path.isdir('./ncbi_fasta'):
        os.chdir('./ncbi_fasta')
    else:
        os.mkdir('./ncbi_fasta')
        os.chdir('./ncbi_fasta')

    wget_command = ['wget', '-O', output_name, ftp_link]
    gunzip_command = ['gunzip', output_name]
    subprocess.run(wget_command)
    subprocess.run(gunzip_command)

    os.chdir('../')


def parse_seqio_fasta(fasta_file, assembly_name):
    """ parses the refseq fasta files for each respective strain, counts
        the number of contigs per SeqIO record, and counts the total
        number of base pairs for each contig from each respective genome
        Args:
            fasta_record_list [lst] : array of SeqIO objects
            log_file [file] : output file needed to write to; opened at main()
        Returns:
            None
    """

    num_of_contigs = 0
    num_of_bp = 0

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        num_of_contigs += 1
        if len(seq_record.seq) >= 1000:
            num_of_bp += len(seq_record.seq)

    with open('UPEC.log', 'a') as log_file:
        log_file.write('There are {} contigs in the {} assembly.\n'.format(\
                                                str(num_of_contigs), str(assembly_name)))

        log_file.write('There are {} base pairs in the {} assembly.\n'.format(\
                                                    str(num_of_bp), str(assembly_name)))

    return None



def main():

    args = arg_parser()

    ftp_files = args.ftp_file

    with open(ftp_files, 'r') as input_ftp:
        for line in input_ftp:
            input_line = line.strip().split(',')

            zipped_fasta_output_name = input_line[0] + '_FASTA.fna.gz'
            zipped_feature_output_name = input_line[0] + '_FEAT.fna.gz'

            wget_gunzip_fasta(zipped_fasta_output_name, input_line[1])
            wget_gunzip_fasta(zipped_feature_output_name, input_line[2])

            fasta_file = "./ncbi_fasta/{}_FASTA.fna".format(input_line[0])

            parse_seqio_fasta(fasta_file, input_line[0])




if __name__ == '__main__':
    main()
