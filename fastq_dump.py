#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0


"""


import argparse
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-s', '--sra_file', help='sra accession id')

    return parser.parse_args()


def prefetch_fastq_decomp(sra_accesion, fastq_dir):
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

    subprocess.run(['prefetch', sra_accesion])

    fq_command = "fastq-dump -I --split-files {} -O {}".format(sra_accesion, fastq_dir)

    subprocess.run(fq_command, shell=True)



def main():
    """ running fastq dump script """
    
    args = arg_parser()

    assembly_name = args.assembly_name
    sra_accesion = args.sra_file

    fastq_dir = assembly_name + '_sra'

    with open(sra_accesion, 'r') as sra_input:
        for line in sra_input:
            sra_list = line.strip().split(',')
            fastq_dir = sra_list[0] + '_sra'
            prefetch_fastq_decomp(sra_list[1], fastq_dir)

if __name__ == '__main__':
    main()
