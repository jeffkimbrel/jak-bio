#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse

import jak_utils
from jakomics import colors, utilities
jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Give it a multifasta input (-f) and a length (-l), it will return a fasta list of sequences that length or longer')

parser.add_argument('-f', '--fasta',
                    help="Fasta input",
                    required=True)

parser.add_argument('-l', '--length',
                    default=1000,
                    help="Length",
                    type=int)

args = parser.parse_args()

# MAIN ########################################################################

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    if len(seq_record.seq) >= args.length:
        print(seq_record.format("fasta").rstrip())
