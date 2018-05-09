#!/usr/bin/python

# Give it a multifasta input and a length, it will return a fasta list of sequences that length or longer

import sys
from Bio import SeqIO

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    if len(seq_record.seq) >= int(sys.argv[2]):
        #print seq_record.format("fasta") # with wrapping
        print(">" + str(seq_record.id) + "\n" + str(seq_record.seq))