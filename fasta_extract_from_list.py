import pandas as pd
from Bio import SeqIO
import argparse
import sys
import os
from jakomics import colors, utilities

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Give it a multifasta input (-f) and a list.txt (-l) with identifiers, and it returns the sequences')

parser.add_argument('--in_dir',
                    help="Directory with fasta files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('-l', '--list',
                    help="List of headers to include",
                    required=True)

parser.add_argument('-o', '--out',
                    help="Fasta out",
                    required=True)

parser.add_argument('-c', '--column',
                    help="0-based column position",
                    type=int,
                    default=0,
                    required=False)

args = parser.parse_args()

# FUNCTIONS ###################################################################

def main(fasta, list, column):
    global found

    original_sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    df = pd.read_csv(list, sep="\t", header=None)
    wanted_sequences = df.iloc[:, column]
    
    for id in wanted_sequences:
        if id in original_sequences:
            found += 1
            collected_sequences.append(original_sequences[id])
        else:
            no.append(id)

# MAIN ########################################################################

if __name__ == "__main__":
    jak_utils.header()


    collected_sequences = []
    found = 0
    no = []

    fasta_files = utilities.get_files(args.files, args.in_dir, ['faa', 'fa', 'ffn', 'fasta'])

    for fasta_file in fasta_files:
        main(fasta_file.file_path, args.list, args.column)

    print(f'{colors.bcolors.GREEN}Found = {found}{colors.bcolors.END}, {colors.bcolors.RED}NOT Found = {len(no)}{colors.bcolors.END}', end="\r", file=sys.stderr)

    print(f'\nWriting results to {args.out}')
    SeqIO.write(collected_sequences, args.out, "fasta")

    if len(no) > 0:
        print(f'Done {colors.bcolors.RED}(with errors){colors.bcolors.END}')
    else:
        print(f'{colors.bcolors.GREEN}Done!{colors.bcolors.END}')

    for id in no:
        print(f'{colors.bcolors.RED}Not Found: {id}{colors.bcolors.END}')

    
