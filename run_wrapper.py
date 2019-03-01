#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

"""
import os
import argparse
import subprocess


CURRENT_DIR = os.getcwd()

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


    ftp_path = CURRENT_DIR + '/' + ftp_files
    sra_path = CURRENT_DIR + '/' + sra_file
    ftp_to_path = CURRENT_DIR + '/Timothy_Baker_tests/' + ftp_files
    sra_to_path = CURRENT_DIR + '/Timothy_Baker_tests/' + sra_file

    # parse_fasta_path = CURRENT_DIR + '/parse_fasta.py'
    # fastq_dump_path = CURRENT_DIR + '/fastq_dump.py'
    # prokka_path = CURRENT_DIR + '/prokka.py'
    # tophat2_path = CURRENT_DIR + '/tophat2.py'
    # cufflinks_path = CURRENT_DIR + '/cufflinks.py'
    #
    # parse_fasta_to_path = CURRENT_DIR + '/Timothy_Baker_tests/parse_fasta.py'
    # fastq_dump_to_path = CURRENT_DIR + '/Timothy_Baker_tests/fastq_dump.py'
    # prokka_to_path = CURRENT_DIR + '/Timothy_Baker_tests/prokka.py'
    # tophat2_to_path = CURRENT_DIR + '/Timothy_Baker_tests/tophat2.py'
    # cufflinks_to_path = CURRENT_DIR + '/Timothy_Baker_tests/cufflinks.py'
    #
    os.rename(ftp_path, ftp_to_path)
    os.rename(sra_path, sra_to_path)
    # os.rename(parse_fasta_path, parse_fasta_to_path)
    # os.rename(fastq_dump_path, fastq_dump_to_path)
    # os.rename(prokka_path, prokka_to_path)
    # os.rename(tophat2_path, tophat2_to_path)
    # os.rename(cufflinks_path, cufflinks_to_path)


    with open(ftp_files, 'r') as ftp_input:
        assembly_name_list = [line.strip().split(',')[0] for line in ftp_input]
        nargs_assembly = '\t'.join(assembly_name_list)

    parse_fasta_cmd = "python3 ../parse_fasta.py -f {}".format(ftp_files)
    subprocess.run(parse_fasta_cmd, shell=True)

    # prokka_cmd = "python3 ../prokka.py -n {}".format(nargs_assembly)
    # subprocess.run(prokka_cmd, shell=True)

    # fastq_dump_cmd = "python3 ../fastq_dump.py -s {}".format(sra_file)
    # subprocess.run(fastq_dump_cmd, shell=True)

    tophat2_cmd = "python3 ../tophat2.py -s {} -t {}".format(sra_file, threads)
    subprocess.run(tophat2_cmd, shell=True)

    cufflinks_cmd = "python3 ../cufflinks.py -n {} -t {}".format(nargs_assembly, threads)
    subprocess.run(cufflinks_cmd, shell=True)







if __name__ == '__main__':
    main()
