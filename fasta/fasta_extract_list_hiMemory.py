#!/usr/bin/python

# Give it a multifasta input and a list.txt with identifiers, and it returns the sequences

# this version goes faster, but uses a lot memory by loading the whole fasta file into memory

import sys
from Bio import SeqIO

## Read in list
lines = [line.strip() for line in open(sys.argv[2])]

handle = open(sys.argv[1], "rU")
record_dict = SeqIO.index("transcripts/transcripts.cd-hit.50aa.faa", "fasta")
handle.close()


for line in list(lines):
    if line in record_dict:
        print(">" + str(record_dict[line].id) + "\n" + str(record_dict[line].seq))