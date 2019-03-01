#!/usr/bin/env python3
"""
@author: Timothy Baker
@version: 1.0.0
cufflinks.py



"""



 #   _____ _    _ ______ ______ _      _____ _   _ _  __ _____
 #  / ____| |  | |  ____|  ____| |    |_   _| \ | | |/ // ____|
 # | |    | |  | | |__  | |__  | |      | | |  \| | ' /| (___
 # | |    | |  | |  __| |  __| | |      | | | . ` |  <  \___ \
 # | |____| |__| | |    | |    | |____ _| |_| |\  | . \ ____) |
 #  \_____|\____/|_|    |_|    |______|_____|_| \_|_|\_\_____/
 #


import argparse
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-n', '--nargs', nargs='+')
    parser.add_argument('-t', '--threads', help='path to fasta file.')

    return parser.parse_args()


def run_cufflinks(gff_file, output_dir, threads, bam_file):
    """ runs the cufflinks software for assembly transcripts
        Args:
            gff_file (str) : path to gff_file
            output_dir (str) : output_dir name
            threads (int) : number of threads to run
            bam_file (str) : path to bam_file
        Returns:
            None
     """

    cufflink_command = "cufflinks -p {} -G {} -o {} {}".format(threads, \
                                                                gff_file, \
                                                                output_dir, \
                                                                bam_file)
    subprocess.run(cufflink_command, shell=True)



def run_cuffmerge(threads, assembly_file):
    """ merges all the assemblies together """

    cuffmerge_command = "cuffmerge -p {} -o {} {}".format(threads, 'merged_ecoli', assembly_file)

    subprocess.run(cuffmerge_command, shell=True)


def run_cuffnorm(threads, assembly_nargs):
    """ runs the cuffnorm software with the merged transcripts from the
        different genomes
        Args:
            threads (int) : number of threads to run
            assembly_nargs (lst) : list of assembly names
    """

    concat_bams = build_bam_list(assembly_nargs)

    merged_gtf = './merged_ecoli/merged.gtf'

    cuffnorm_command = "cuffnorm -o diff_results -p {} {} {}".format(threads, \
                                                        merged_gtf, concat_bams)

    subprocess.run(cuffnorm_command, shell=True)


def build_cuff_run(assembly_name_list, threads):
    """ builds the entire run, and sets the path names
        Args:
            assembly_name_list (lst) : list of assembly names from input
            threads (int) : number of threads to run
        Returns:
            None
     """

    for assembly_name in assembly_name_list:

        gff_file = './' + assembly_name + '_prokout/' + assembly_name + '_index.gff'
        output_dir = assembly_name + '_cuff'
        bam_file = assembly_name + '_tophat/' + 'accepted_hits.bam'

        run_cufflinks(gff_file, output_dir, threads, bam_file)

        with open('assemblies.txt', 'a') as assemble_file:
            assemble_file.write("./{}_cuff/transcripts.gtf\n".format(assembly_name))



def build_bam_list(assembly_nargs):
    """ creates a list of bam paths to concatenate together for easy placement
        into command line
        Args:
            assembly_nargs (lst) : list of strain names
        Returns:
            concatenated bam list by tab
    """
    bam_list = []

    for assemble_name in assembly_nargs:
        bam_file = assemble_name + '_tophat/' + 'accepted_hits.bam'
        bam_list.append(bam_file)

    return '\t'.join(bam_list)

ls -
def main():
    """ runs main script """

    # sets parser
    args = arg_parser()

    # sets arguments
    assembly_nargs = list(args.nargs)
    threads = args.threads

    # runs the build with the nargs
    # nargs need to be in list format
    build_cuff_run(assembly_nargs, threads)

    # runs the cuffmerge
    run_cuffmerge(threads, 'assemblies.txt')

    # runs cuffnorm
    run_cuffnorm(threads, assembly_nargs)

if __name__ == '__main__':
    main()
