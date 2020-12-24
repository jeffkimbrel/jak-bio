#!/usr/bin/env python

from Bio import SeqIO
import argparse
import os
import pandas as pd

from jakomics import colors
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Give a fasta file and a grouping file, and it will split into smaller fasta files with group as file name')

parser.add_argument('--fasta',
                    help="Fasta file",
                    required=True)

parser.add_argument('-g', '--group_file',
                    help="Paths to tsv file with sequence group info. Column 1 is fasta header, column 2 is the group",
                    required=True)

parser.add_argument('--out_dir',
                    help="Directory to write individual fasta files to",
                    required=True)

args = parser.parse_args()

# FUNCTIONS ###################################################################


def get_fasta():
    seqs = {}
    for seq_record in SeqIO.parse(args.fasta, "fasta"):
        seqs[seq_record.id] = seq_record
    return seqs


def parse_group(header_col=0, group_col=1):
    groups = {}

    df = pd.read_csv(args.group_file, sep="\t", header=None)

    for index, row in df.iterrows():
        if row[group_col] in groups:
            groups[row[group_col]].append(row[header_col])
        else:
            groups[row[group_col]] = [row[header_col]]

    return groups


# MAIN ########################################################################
if __name__ == "__main__":
    jak_utils.header()

    args.out_dir = os.path.abspath(args.out_dir) + '/'

    if not os.path.exists(args.out_dir):
        print(f"{colors.bcolors.BLUE}Creating directory {args.out_dir}{colors.bcolors.END}")
        os.makedirs(args.out_dir)

    original_seqs = get_fasta()
    groups = parse_group()

    for group in groups:

        group_seqs = []

        for header in groups[group]:
            if header in original_seqs:
                group_seqs.append(original_seqs[header])
            else:
                raise Exception(
                    f"{colors.bcolors.RED}ERROR: {header} is not found in {args.fasta}{colors.bcolors.END}")

        SeqIO.write(group_seqs, os.path.join(args.out_dir, group + '.fa'), "fasta")
