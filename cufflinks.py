#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0
cufflinks.py



"""

import argparse
import subprocess

 #   _____ _    _ ______ ______ _      _____ _   _ _  __ _____
 #  / ____| |  | |  ____|  ____| |    |_   _| \ | | |/ // ____|
 # | |    | |  | | |__  | |__  | |      | | |  \| | ' /| (___
 # | |    | |  | |  __| |  __| | |      | | | . ` |  <  \___ \
 # | |____| |__| | |    | |    | |____ _| |_| |\  | . \ ____) |
 #  \_____|\____/|_|    |_|    |______|_____|_| \_|_|\_\_____/
 #



def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-n', '--nargs', nargs='+')
    parser.add_argument('-t', '--threads', help='path to fasta file.')

    return parser.parse_args()


def run_cufflinks(gff_file, output_dir, threads, bam_file):

    cufflink_command = "cufflinks -p {} -G {} -o {} {}".format(threads, \
                                                                gff_file, \
                                                                output_dir, \
                                                                bam_file)
    subprocess.run(cufflink_command, shell=True)



def run_cuffmerge(threads, assembly_file):

    cuffmerge_command = "cuffmerge -p {} -o {} {}".format(threads, 'merged_ecoli', assembly_file)
    subprocess.run(cuffmerge_command, shell=True)



# def run_cuffnorm(threads, merged_gtf, bam_list):
#     cuffnorm_command = "cuffnorm -o diff_results -p {} {} {} {} {} {}".format(threads, merged_gtf)
#
#     subprocess.run(cuffnorm_command, shell=True)

def build_cuff_run(assembly_name_list, threads):
    for assembly_name in assembly_name_list:

        gff_file = './' + assembly_name + '_prokout/' + assembly_name + '_index.gff'
        output_dir = assembly_name + '_cuff'
        bam_file = assembly_name + '_tophat/' + 'accepted_hits.bam'
        run_cufflinks(gff_file, output_dir, threads, bam_file)

        with open('assemblies.txt', 'a') as assemble_file:
            assemble_file.write("./{}_cuff/transcripts.gtf\n".format(assembly_name))


def main():

    args = arg_parser()

    nargs = args.nargs
    threads = args.threads

    build_cuff_run(nargs, threads)

    run_cuffmerge(threads, 'assemblies.txt')

if __name__ == '__main__':
    main()
