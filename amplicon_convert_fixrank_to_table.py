#!/usr/bin/env python

import pandas as pd
import argparse
import sys
from jakomics.taxa import RDP

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Convert an RDP fixrank table to a taxa table')

parser.add_argument('-f', '--file',
                    help="RDP fixrank file",
                    required=True)

parser.add_argument('-t', '--threshold',
                    help="RDP confidence threshold",
                    type=float,
                    default=50)

args = parser.parse_args()

if args.threshold > 1:
    args.threshold /= 100

# MAIN ########################################################################


if __name__ == "__main__":
    jak_utils.header()

    column_names = ['ASV', 'STRAND', 'D', 'D_p', 'P', 'P_p',
                    'C', 'C_p', 'O', 'O_p', 'F', 'F_p', 'G', 'G_p']

    df = pd.read_csv(args.file,
                     skiprows=7,
                     sep=";",
                     header=None,
                     names=column_names,
                     index_col=None)

    print("", "domain", "phylum", "class", "order", "family", "genus", sep="\t")

    unclassified_counts = {'domain': 0, 'phylum': 0,
                           'class': 0, 'order': 0, 'family': 0, 'genus': 0}

    for id, row in df.iterrows():
        asv = RDP(row, args.threshold)
        asv.view()
        unclassified_counts = asv.count_unclassified(unclassified_counts)

    print(unclassified_counts, file=sys.stderr)
