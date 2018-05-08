import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='xxx')

parser.add_argument('-f', '--fasta',
    help="fasta file",
    required=True)

parser.add_argument('-p', '--prefix',
    help="Prefix",
    required=True)

parser.add_argument('--sort', '-s',
    action = 'store_true',
    help = 'Sort from largest to smallest before re-numbering' )

args = parser.parse_args()

seqs = {}
lengths = {}

counter = 1

for record in SeqIO.parse(args.fasta, "fasta"):
    seqs[str(record.id)] = str(record.seq).upper()
    lengths[str(record.id)] = len(str(record.seq).upper())

oldHeaders = lengths.keys()

if args.sort:
    oldHeaders = sorted(lengths, key = lengths.get)
    oldHeaders.reverse()

for head in oldHeaders:
    print(">" + args.prefix + "_" + str(counter), seqs[head], sep = "\n")
    counter += 1
