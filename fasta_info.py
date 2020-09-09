#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
from Bio.SeqUtils import GC

import jak_utils
from jakomics import colors, utilities
jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Given a (list of) fasta file(s), returns the total or subsequence nt/aa length and GC%.')

parser.add_argument('--in_dir',
                    help="Directory with fasta files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--total',
                    '-t',
                    action='store_true',
                    help='Print total of all subsequences')

parser.add_argument('--sequence',
                    '-s',
                    action='store_true',
                    help='Print total for each subsequence')

args = parser.parse_args()

if args.total is False and args.sequence is False:
    print(f"{colors.bcolors.RED}ERROR: at least one of --total and --sequence required{colors.bcolors.END}", file=sys.stderr)
    sys.exit()

# FUNCTIONS ###################################################################


def getInfo(seq_record):
    id = str(seq_record.id)
    count = len(str(seq_record.seq))
    return(id, count)

# ANALYZE #####################################################################


fasta_files = utilities.get_files(args.files, args.in_dir, ['faa', 'fa', 'ffn', 'fasta'])

for fasta_file in fasta_files:
    total = 0
    total_gc = 0
    for seq_record in SeqIO.parse(fasta_file.file_path, "fasta"):

        id, count = getInfo(seq_record)
        total += count
        seq_gc = sum(seq_record.seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])
        total_gc += seq_gc

        if args.sequence is True:
            print(fasta_file.name, id, count, seq_gc / count, sep="\t")

    if args.sequence is False:
        print(fasta_file.name, total, total_gc / total, sep="\t")

if args.total is True and args.sequence is True:
    print(f"{colors.bcolors.YELLOW}WARNING: if both '-s' and '-t' are given, only the sequence lengths will be returned.{colors.bcolors.END}", file=sys.stderr)
