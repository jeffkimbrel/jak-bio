import pandas as pd
from Bio import SeqIO
import argparse
import sys
import os
from jakomics import colors

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Give it a multifasta input (-f) and a list.txt (-l) with identifiers, and it returns the sequences')

parser.add_argument('-f', '--fasta',
                    help="Fasta input",
                    required=True)

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


def main(fasta, list, column, out_file):
    print(f'Reading and processing {fasta}')
    original_sequences = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    collected_sequences = []

    print(f'Reading and processing {list} column {column}')
    df = pd.read_csv(list, sep="\t", header=None)
    wanted_sequences = df.iloc[:, column]

    yes = 0
    no = []

    for id in wanted_sequences:
        if id in original_sequences:
            yes += 1
            collected_sequences.append(original_sequences[id])
        else:
            no.append(id)

        print(f'IDs in {os.path.basename(fasta)}: {colors.bcolors.GREEN}Found = {yes}{colors.bcolors.END}, {colors.bcolors.RED}NOT Found = {len(no)}{colors.bcolors.END}', end="\r", file=sys.stderr)

    print(f'\nWriting results to {out_file}')
    SeqIO.write(collected_sequences, out_file, "fasta")

    if len(no) > 0:
        print(f'Done {colors.bcolors.RED}(with errors){colors.bcolors.END}')
    else:
        print(f'{colors.bcolors.GREEN}Done!{colors.bcolors.END}')

    for id in no:
        print(f'{colors.bcolors.RED}Not Found: {id}{colors.bcolors.END}')


# MAIN ########################################################################

#
# lines = [line.strip() for line in open(args.list)]
#
# handle = open(args.fasta, "rU")
#
# for seq_record in SeqIO.parse(handle, "fasta"):
#     for line in list(lines):
#         if line == seq_record.id:
#             print(">" + str(seq_record.id) + "\n" + str(seq_record.seq))
#             lines.remove(line)


if __name__ == "__main__":
    jak_utils.header()

    main(args.fasta, args.list, args.column, args.out)
