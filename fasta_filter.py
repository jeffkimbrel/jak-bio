#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import re
import os
from tqdm import tqdm

from jakomics import colors, utilities
import jak_utils
jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Filters fasta files, writing to a new file appended with ".ff"')

parser.add_argument('--in_dir',
                    help="Directory with fasta files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--min_length',
                    default=1,
                    help="Length",
                    type=int)

parser.add_argument('--remove_trailing_asterisks',
                    action='store_true',
                    help='Remove stop codon asterisks')

args = parser.parse_args()

# FUNCTIONS ###################################################################


def remove_trailing_asterisks(seqs):
    asterisks_removed = []

    for seq in seqs:
        if seq.seq.endswith("*"):
            seq.seq = Seq(re.sub(r"\*$", "", str(seq.seq)))

        asterisks_removed.append(seq)

    return asterisks_removed


def length_check(seqs):
    length_checked = []

    for seq in seqs:
        if len(str(seq.seq)) >= args.min_length:
            length_checked.append(seq)

    return length_checked


def write_fasta(seqs, ff):
    base_dir = os.path.dirname(ff.file_path)
    # for seq in seqs:
    #     print(seq.id, seq.seq)
    SeqIO.write(seqs, os.path.join(base_dir, ff.new_name), "fasta")

# MAIN ########################################################################


fasta_files = utilities.get_files(args.files, args.in_dir, ['faa', 'fa', 'ffn', 'fasta'])

pbar = tqdm(total=len(fasta_files), desc="Finished", unit=" fasta files")

for fasta_file in fasta_files:
    fasta_file.new_name = f'{fasta_file.short_name}.ff{fasta_file.suffix}'

    seqs = []
    for seq_record in SeqIO.parse(fasta_file.file_path, "fasta"):
        seqs.append(seq_record)

    if args.remove_trailing_asterisks:
        seqs = remove_trailing_asterisks(seqs)

    seqs = length_check(seqs)

    write_fasta(seqs, fasta_file)
    pbar.update(1)
