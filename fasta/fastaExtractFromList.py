#!/usr/bin/python

import sys
from Bio import SeqIO
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'Give it a multifasta input and a list.txt with identifiers, and it returns the sequences')

parser.add_argument('-f', '--fasta',
    help = "Fasta input",
    required = True)

parser.add_argument('-l', '--list',
    help = "List of headers to include",
    required = True)

args = parser.parse_args()

## Read in list
lines = [line.strip() for line in open(args.list)]

handle = open(args.fasta, "rU")

for seq_record in SeqIO.parse(handle, "fasta"):

    for line in list(lines):
        if line == seq_record.id:
            print(">" + str(seq_record.id) + "\n" + str(seq_record.seq))
            lines.remove(line)
