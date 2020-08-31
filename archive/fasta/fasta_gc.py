import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

from jakomics import utilities, colors

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with table files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual table files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--total',
                    '-t',
                    action='store_true',
                    help='Calculate total for all subsequences in a file')

parser.add_argument('--sequence',
                    '-s',
                    action='store_true',
                    help='Calculate individually for all subsequences in a file')

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == '__main__':

    if args.total is False:
        if args.sequence is False:
            print(
                f"{colors.bcolors.RED}You should probably use the -s or -t flags to see some results...{colors.bcolors.END}")

            sys.exit()

    file_list = utilities.get_files(args.files, args.in_dir, ["fa", "fna", "fasta"])
    for file in file_list:
        total_length = 0
        total_gc = 0
        for seq_record in SeqIO.parse(file.file_path, "fasta"):
            seq = seq_record.seq
            seq_length = len(seq)
            seq_gc = sum(seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])

            if args.sequence:
                print(f'{file.short_name}\t{seq_length}\t{100*seq_gc / seq_length}')
            else:
                total_length += seq_length
                total_gc += seq_gc

        if args.total is True:
            print(f'{file.short_name}\t{total_length}\t{100*total_gc / total_length}')
