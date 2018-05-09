#!/usr/bin/python

# Give it a multifasta input and a list.txt with identifiers, and it returns the sequences

# this version goes slower, but uses less memory, by not loading the whole fasta file into memory

import sys
from Bio import SeqIO

## Read in list
lines = [line.strip() for line in open(sys.argv[2])]

handle = open(sys.argv[1], "rU")

for seq_record in SeqIO.parse(handle, "fasta"):
    
    for line in list(lines):
        if line == seq_record.id:             
            print(">" + str(seq_record.id) + "\n" + str(seq_record.seq))
            lines.remove(line)