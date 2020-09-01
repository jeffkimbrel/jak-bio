#!/usr/bin/env python

import pandas as pd
import argparse
import sys
from Bio import SeqIO

import jak_utils
from jakomics import colors
jak_utils.header()

# OPTIONS #####################################################################
parser = argparse.ArgumentParser(description='TPM transform columns in a table')

parser.add_argument('-l', '--lengths',
                    help="text file with 'feature' and 'length' columns",
                    required=False)

parser.add_argument('-f', '--fasta',
                    help="fasta file to calculate the lengths from",
                    required=False)

parser.add_argument('-t', '--table',
                    help="Table with raw counts",
                    required=True)

parser.add_argument('-o', '--out',
                    help="Filename to write TPM table to",
                    required=True)

args = parser.parse_args()

if args.lengths is None and args.fasta is None:
    print(f"\n{colors.bcolors.RED}ERROR: at least one of --lengths and --fasta required{colors.bcolors.END}")
    sys.exit()

# FUNCTIONS ###################################################################


def get_lengths():
    if args.lengths is not None:
        lengths = pd.read_csv(args.lengths, sep="\t")
        lengths = lengths.set_index('feature')
        lengths = lengths['length']
    else:
        lengths = pd.DataFrame(columns=['feature', 'length'])
        for seq_record in SeqIO.parse(args.fasta, "fasta"):
            id = str(seq_record.id)
            count = len(str(seq_record.seq))
            lengths = lengths.append({'feature': id, 'length': count}, ignore_index=True)
        lengths = lengths.set_index('feature')
    print(f'Lengths file has {len(lengths)} features')
    return lengths


def prepare_table():
    table = pd.read_csv(args.table, sep="\t")
    table = table.set_index(list(table.columns[[0]]))
    table = table.fillna(0)
    print(f'Counts table has {len(table)} features')
    return table


def merge_tables(table, lengths):
    df = pd.concat([table, lengths], join='inner', axis=1)
    df.index.name = table.index.name
    print(f'Final table merged {len(df)} common features')
    return df


def calc_rpk(df):
    rpk = df.div(df.length, axis=0)
    rpk = rpk.drop(columns=['length'])
    return rpk


def calc_tpm(rpk, per_million):
    tpm = rpk.div(per_million, axis=1)
    tpm.index.name = rpk.index.name
    return tpm

# MAIN ########################################################################


if __name__ == "__main__":
    lengths = get_lengths()
    table = prepare_table()
    df = merge_tables(table, lengths)
    rpk = calc_rpk(df)
    per_million = rpk.sum(axis=0) / 1000000
    tpm = calc_tpm(rpk, per_million)

    if tpm.index.name == "Unnamed: 0":
        tpm.index.name = "TPM"
    else:
        tpm.index.name = tpm.index.name + "_TPM"

    tpm.to_csv(args.out, sep="\t")

    # print("\nPreview:")
    # print(tpm)

    print("\nTPM adjusted values saved to " + args.out)
