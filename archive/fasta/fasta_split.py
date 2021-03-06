#!/usr/bin/python

import sys
from Bio import SeqIO
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Give it a multifasta input and an integer "x", and it will make as many smaller files necessary with at most that x sequences.')

parser.add_argument('-f', '--fasta',
                    help="Fasta input",
                    required=True)

parser.add_argument('-s', '--sequences',
                    help="Maximum number of sequences per file",
                    required=True)

parser.add_argument('-o', '--out',
                    default="out",
                    help="Output file names (appended with an integer and \".fa\")")

args = parser.parse_args()
args.sequences = int(args.sequences)

## MAIN ########################################################################

seqCounter = 0
fileCounter = 1
fileTarget = open(args.out+""+str(fileCounter)+".fa", 'w')

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    if seqCounter == args.sequences:
        fileTarget.close()
        fileCounter += 1
        fileTarget = open(args.out+""+str(fileCounter)+".fa", 'w')
        seqCounter = 0

    seqCounter += 1

    fileTarget.write(">"+seq_record.id+"\n"+str(seq_record.seq)+"\n")
    # seq = str(seq_record.seq)
    # seq = seq.replace("-", "")
    # fileTarget.write(">"+args.out+""+str(fileCounter)+"\n"+seq+"\n")
    # print(seq_record.id, args.out+""+str(fileCounter), sep="\t")

fileTarget.close()
