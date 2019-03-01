#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0

"""
import os
import logging
import argparse
import subprocess


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



CURRENT_DIR = os.getcwd()

if os.path.isdir('./Timothy_Baker'):
    os.chdir('./Timothy_Baker')
else:
    os.mkdir('./Timothy_Baker')
    os.chdir('./Timothy_Baker')


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
    """ runs main script """

    LOGGER.info("Beginning to run wrapper.")

    args = arg_parser()
    threads = args.threads
    ftp_files = args.ftp_links
    sra_file = args.sra_file

    LOGGER.info("The arguments input are %s", str(args))

    ftp_path = CURRENT_DIR + '/' + ftp_files
    sra_path = CURRENT_DIR + '/' + sra_file
    ftp_to_path = CURRENT_DIR + '/Timothy_Baker/' + ftp_files
    sra_to_path = CURRENT_DIR + '/Timothy_Baker/' + sra_file

    LOGGER.info("Moving input files into new directory.")
    os.rename(ftp_path, ftp_to_path)
    os.rename(sra_path, sra_to_path)

    with open(ftp_files, 'r') as ftp_input:
        assembly_name_list = [line.strip().split(',')[0] for line in ftp_input]
        nargs_assembly = '\t'.join(assembly_name_list)
        LOGGER.info("The assembly names are %s", str(nargs_assembly))

    LOGGER.info("Downloading FASTA and Feature files from NCBI.")
    parse_fasta_cmd = "python3 ../scripts/parse_fasta.py -f {}".format(ftp_files)
    subprocess.run(parse_fasta_cmd, shell=True)

    LOGGER.info("Running prokkaself.")
    prokka_cmd = "python3 ../scripts/prokka.py -n {}".format(nargs_assembly)
    subprocess.run(prokka_cmd, shell=True)

    LOGGER.info("Downloading SRA files and decompressing to FASTQ")
    # fastq_dump_cmd = "python3 ../scripts/fastq_dump.py -s {}".format(sra_file)
    # subprocess.run(fastq_dump_cmd, shell=True)

    # LOGGER.info("Building bowtie2 index and running Tophat2 alignment.")
    # tophat2_cmd = "python3 ../scripts/tophat2.py -s {} -t {}".format(sra_file, threads)
    # subprocess.run(tophat2_cmd, shell=True)

    LOGGER.info("Assembling transcripts, merging, and normalization.")
    cufflinks_cmd = "python3 ../scripts/cufflinks.py -n {} -t {}".format(nargs_assembly, threads)
    subprocess.run(cufflinks_cmd, shell=True)

    LOGGER.info("Pipeline complete.")


if __name__ == '__main__':
    main()
