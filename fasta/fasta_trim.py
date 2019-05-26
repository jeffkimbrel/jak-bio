from Bio import SeqIO
import sys
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='takes a multifasta file and trims from the beginning and end by a specified amount')

parser.add_argument('-f', '--fasta',
    help="fasta file",
    required=True)

parser.add_argument('-b', '--beginning',
    default=0,
    help="Trim from the beginning (5')" )

parser.add_argument('-e', '--end',
    default=0,
    help="Trim from the end (3')" )

args = parser.parse_args()

## MAIN ########################################################################

handle = open(args.fasta, "rU")
for record in SeqIO.parse(handle, "fasta"):
    seq = record.seq

    trimmedSeq = ""

    if args.end == 0:
        trimmedSeq = seq[int(args.beginning):]
    else:
        trimmedSeq = seq[int(args.beginning):(0-int(args.end))]

    #trimmedSeq = seq[0:]

    #print(seq,trimmedSeq,sep="\n")
    print(">"+str(record.id)+"\n"+trimmedSeq)
