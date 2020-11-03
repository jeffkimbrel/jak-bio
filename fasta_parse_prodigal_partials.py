import os
import sys
import argparse
from Bio import SeqIO

from jakomics import utilities, colors
import jak_utils

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

parser.add_argument('-p', '--partial',
                    help="Partial code(s) to keep",
                    default=['00'],
                    nargs='*',
                    type=str)

parser.add_argument('--strip',
                    '-s',
                    action='store_true',
                    help='Strip gene header to just ID')

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == '__main__':

    jak_utils.header()

    file_list = utilities.get_files(args.files, args.in_dir, ["fa", "fna", "fasta", "faa"])

    for file in file_list:
        file.out_path = file.name + ".partial_parsed" + file.suffix
        with open(file.out_path, 'w') as out_f:
            for record in SeqIO.parse(file.file_path, 'fasta'):
                for partial in args.partial:
                    if 'partial=' + str(partial) in record.description:
                        if args.strip:
                            record.description = record.id

                        SeqIO.write(record, out_f, 'fasta')
