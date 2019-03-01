#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

"""
import os
import argparse
import subprocess


CURRENT_DIR = os.getcwd

if os.path.isdir('./Timothy_Baker_tests'):
    os.chdir('./Timothy_Baker_tests')
else:
    os.mkdir('./Timothy_Baker_tests')
    os.chdir('./Timothy_Baker_tests')


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-p', '--ftp_links', help='path to ftp links input.')
    parser.add_argument('-s', '--sra_file', help='text file of sra ids')
    parser.add_argument('-t', '--threads', help='number of threads to use')
    return parser.parse_args()



def main():

    args = arg_parser()

    threads = args.threads
    ftp_files = args.ftp_links
    sra_file = args.sra_file

    os.rename(CURRENT_DIR + ftp_files, './Timothy_Baker_tests/' + ftp_files)
    os.rename(CURRENT_DIR + sra_file, './Timothy_Baker_tests/' + sra_file)

    with open(ftp_files, 'r') as ftp_input:
        assembly_name_list = [line.strip().split(',')[0] for line in ftp_input]


    parse_fasta_cmd = "python3 parse_fasta.py -f {}".format(ftp_files)

    subprocess.run(parse_fasta_cmd, shell=True)

    fastq_dump_cmd = "python3 fastq_dump.py -s {}".format(sra_file)

    subprocess.run(fastq_dump_cmd, shell=True)

    prokka_cmd = "python3 prokka.py -a {}".format(assembly_name_list)
    subprocess.run(prokka_cmd, shell=True)

    tophat2_cmd = "python3 tophat2.py -a {} -s {} -t {}".format(assembly_name_list, \
                                                                sra_file, threads)
    subprocess.run(tophat2_cmd, shell=True)

    cufflinks_cmd = "python3 cufflinks.py -a {} -t {}".format(assembly_name_list, threads)
    subprocess.run(cufflinks_cmd, shell=True)







if __name__ == '__main__':
    main()
