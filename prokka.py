#!/usr/bin/env python3
"""
@author: Timothy Baker
@date: 02-25-2019

run_prokka.py

"""

 #  _____  _____   ____  _  ___  __
 # |  __ \|  __ \ / __ \| |/ / |/ /    /\
 # | |__) | |__) | |  | | ' /| ' /    /  \
 # |  ___/|  _  /| |  | |  < |  <    / /\ \
 # | |    | | \ \| |__| | . \| . \  / ____ \
 # |_|    |_|  \_\\____/|_|\_\_|\_\/_/    \_\
 #


import argparse
import shlex
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-a', '--assembly_name', help='base name of the strain file from prokka.')
    parser.add_argument('-f', '--fasta_file', help='path to fasta file.')
    parser.add_argument('-r', '--refseq_file', help='path to refseq feature file.')

    return parser.parse_args()



def run_prokka(fasta_file, output_dir, prefix_name):

    prokka_cmd = "prokka --outdir {} --prefix {} {} --genus Escherichia".format(output_dir, \
                                                                                    prefix_name, \
                                                                                    fasta_file)
    args = shlex.split(prokka_cmd)

    prokka_run = subprocess.Popen(args, stdout=subprocess.PIPE, \
                                stderr=subprocess.PIPE, universal_newlines=True)

    prokka_communicate = prokka_run.communicate()

    return prokka_communicate



def copy_prokka_text(prefix_name, log_file_name):

    copy_cmd = "cat {}.txt >> {}.log".format(prefix_name, log_file_name)

    args = shlex.split(copy_cmd)

    copy_process = subprocess.Popen(args, stdout=subprocess.PIPE, \
                            stderr=subprocess.PIPE, universal_newlines=True)

    copy_process = copy_process.communicate()

    return copy_process


def grep_count(word, input_file):
    """ takes a word, and an input file and returns the count using grep
        Args:
            word (str) : string
            input_file (file) : any input to grep
        Return
            integer object count of one word

     """

    grep_process = subprocess.Popen(["grep", "-c", word, input_file], \
                                                    stdout=subprocess.PIPE, \
                                                    stderr=subprocess.PIPE, \
                                                    universal_newlines=True)

    word_count = grep_process.communicate()

    if word_count:
        return int(word_count[0].strip())

    return None


def cds_trna_difference(prokka_file, refseq_file):

    prokka_dict_count = {'\tCDS\t':0, '\ttRNA\t':0}
    refseq_dict_count = {'\tCDS\t':0, '\ttRNA\t':0}

    for key in prokka_dict_count:
        count = grep_count(key, prokka_file)
        prokka_dict_count[key] = int(count)

    with open(refseq_file, 'r') as refseq:
        for line in refseq:
            new_line = line.strip().split('\t')

            if new_line[0] == 'CDS' and new_line[1] == 'with_protein':
                refseq_dict_count['\tCDS\t'] = int(new_line[5])

            if new_line[0] == 'tRNA':
                refseq_dict_count['\ttRNA\t'] = int(new_line[6])

    result_dict = {}

    for key, value in refseq_dict_count.items():
        result = prokka_dict_count[key] - value
        result_dict[key] = result

    return result_dict.values()


def write_output_tmp(count, name):

    with open('tmp.txt', 'w') as tmp_file:

        if count[0] < 0 and count[1] == 0:
            tmp_file.write("Prokka found {} less CDS and the same amount of tRNA than the Refseq in assembly {}".format(count[0], name))
        elif count[0] > 0 and count[1] == 0:
            tmp_file.write("Prokka found {} additional CDS and the same amount of tRNA than the Refseq in assembly {}".format(count[0], name))

        elif count[0] == 0 and count[1] < 0:
            tmp_file.write("Prokka found the same amount of CDS and {} less tRNA than the Refseq in assembly {}".format(count[1], name))
        elif count[0] == 0 and count[1] > 0:
            tmp_file.write("Prokka found the same amount of CDS and {} additional tRNA than the Refseq in assembly {}".format(count[1], name))

        elif count[0] > 0 and count[1] > 0:
            tmp_file.write("Prokka found {} additional CDS and {} additional tRNA than the Refseq in assembly {}".format(count[0], count[1], name))
        elif count[0] < 0 and count[1] < 0:
            tmp_file.write("Prokka found {} less CDS and {} less tRNA than the Refseq in assembly {}".format(count[0], count[1], name))

        elif count[0] > 0 and count[1] < 0:
            tmp_file.write("Prokka found {} additional CDS and {} less tRNA than the Refseq in assembly {}".format(count[0], count[1], name))
        elif count[0] < 0 and count[1] > 0:
            tmp_file.write("Prokka found {} less CDS and {} additional tRNA than the Refseq in assembly {}".format(count[0], count[1], name))

    return None


def main():

    args = arg_parser()

    assembly_name = args.assembly_name
    fasta_file = args.fasta_file
    refseq_file = args.refseq_file

    prefix_name = assembly_name + '_index'
    output_dir = assembly_name + '_prokout'
    prokka_file = './' + output_dir + '/' + prefix_name + '.gff'


    run_prokka(fasta_file, output_dir, prefix_name)

    count = cds_trna_difference(prokka_file, refseq_file)

    write_output_tmp(count, assembly_name)

    copy_prokka_text('tmp', 'UPEC')

    copy_prokka_text(prokka_file, 'UPEC')


if __name__ == '__main__':
    main()
