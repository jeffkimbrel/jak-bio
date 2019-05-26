#!/usr/bin/python

import sys
from Bio import SeqIO
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = 'Give it a multifasta input (-f) and a length )-l), it will return a fasta list of sequences that length or longer')

parser.add_argument('-f', '--fasta',
    help = "Fasta input",
    required = True)

parser.add_argument('-l', '--length',
    help = "Length",
    type = int,
    required = True)

args = parser.parse_args()

## MAIN ########################################################################

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    if len(seq_record.seq) >= args.length:
        #print seq_record.format("fasta") # with wrapping
        print(">" + str(seq_record.id) + "\n" + str(seq_record.seq))
