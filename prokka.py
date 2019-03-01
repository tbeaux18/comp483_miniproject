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
import subprocess


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='RNA-seq python wrapper for 4 E.Coli genomes.\n'
    )
    parser.add_argument('-n', '--nargs', nargs='+')
    return parser.parse_args()



def run_prokka(fasta_file, prefix_name, output_dir):
    """ runs the prokka software pipeline; UPEC log must be created for
        logging of the command.
        Args:
            fasta_file (str) : path to fasta file
            prefix_name (str) : assembly_name_index
            output_dir (str) : assembly_name_prokout
        Returns:
            None
    """

    prokka_cmd = "prokka --outdir {} --prefix {} {} --genus Escherichia --usegenus".format(output_dir, \
                                                                                    prefix_name, \
                                                                                    fasta_file)
    subprocess.run(prokka_cmd, shell=True)

    with open('UPEC.log', 'a') as output_prokka_command:
        output_prokka_command.write('\n'+str(prokka_cmd))


def copy_prokka_text(prefix_name, log_file_name):
    """ takes a tmp file and copies to log using bash """

    copy_cmd = "cat {}.txt >> {}.log".format(prefix_name, log_file_name)

    subprocess.run(copy_cmd, shell=True)



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
        return word_count[0].strip()

    return None


def cds_trna_difference(prokka_file, refseq_file):
    """ calculates the difference between prokka output and refseq feature
        Args:
            prokka_file (file) : path to file
            refseq_file (file) : path to file
        Returns:
            result_dict.values() [CDS, tRNA]
    """

    prokka_dict_count = {'\tCDS\t':0, '\ttRNA\t':0}
    refseq_dict_count = {'\tCDS\t':0, '\ttRNA\t':0}

    for key in prokka_dict_count:
        count = grep_count(key, prokka_file)
        try:
            prokka_dict_count[key] = int(count)
        except ValueError:
            pass

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
    """ writes to a tmp file so we can append output into the UPEC logself.
        Args:
            count (lst) : lst of counts [CDS, tRNA]
            name (str) : assembly name
        Returns:
            None
    """

    with open('UPEC.log', 'a') as tmp_file:

        if count[0] < 0 and count[1] == 0:
            tmp_file.write("\nProkka found {} less CDS and the same amount of tRNA than the Refseq in assembly {}\n".format(count[0], name))
        elif count[0] > 0 and count[1] == 0:
            tmp_file.write("\nProkka found {} additional CDS and the same amount of tRNA than the Refseq in assembly {}\n".format(count[0], name))

        elif count[0] == 0 and count[1] < 0:
            tmp_file.write("\nProkka found the same amount of CDS and {} less tRNA than the Refseq in assembly {}\n".format(count[1], name))
        elif count[0] == 0 and count[1] > 0:
            tmp_file.write("\nProkka found the same amount of CDS and {} additional tRNA than the Refseq in assembly {}\n".format(count[1], name))

        elif count[0] > 0 and count[1] > 0:
            tmp_file.write("\nProkka found {} additional CDS and {} additional tRNA than the Refseq in assembly {}\n".format(count[0], count[1], name))
        elif count[0] < 0 and count[1] < 0:
            tmp_file.write("\nProkka found {} less CDS and {} less tRNA than the Refseq in assembly {}\n".format(count[0], count[1], name))

        elif count[0] > 0 and count[1] < 0:
            tmp_file.write("\nProkka found {} additional CDS and {} less tRNA than the Refseq in assembly {}\n".format(count[0], count[1], name))
        elif count[0] < 0 and count[1] > 0:
            tmp_file.write("\nProkka found {} less CDS and {} additional tRNA than the Refseq in assembly {}\n".format(count[0], count[1], name))

    return None


def build_prokka_run(assembly_name_list):
    """ assembles and sets all the file names and bgins to run the module.
        Args:
            assembly_name_list (lst) : list of assembly names
        Returns:
            None
    """

    for assembly_name in assembly_name_list:

        fasta_file = './ncbi_fasta/' + assembly_name + '_FASTA.fna'
        prefix_name = assembly_name + '_index'
        output_dir = assembly_name + '_prokout'
        prokka_file = './' + output_dir + '/' + prefix_name + '.gff'
        prokka_text = './' + output_dir + '/' + prefix_name + '.txt'
        refseq_file = './ncbi_fasta/' + assembly_name + '_FEAT.fna'

        run_prokka(fasta_file, prefix_name, output_dir)

        count = cds_trna_difference(prokka_file, refseq_file)

        write_output_tmp(list(count), assembly_name)

        with open(prokka_text, 'r') as input_file:
            with open('UPEC.log', 'a') as output_file:
                for line in input_file:
                    output_file.write('\n' + line + '\n')

def main():
    """ runs mains script """

    # sets parser
    args = arg_parser()

    # assembly must be in list format
    assembly_nargs = list(args.nargs)

    # runs the build; no error handling available
    build_prokka_run(assembly_nargs)


if __name__ == '__main__':
    main()
