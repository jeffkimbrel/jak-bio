#!/usr/bin/env python

import argparse
import pandas as pd
from tqdm import tqdm

from jakomics.fastq import FASTQ
from jakomics import colors

import jak_utils


# OPTIONS #####################################################################
parser = argparse.ArgumentParser(
    description='test')

parser.add_argument('-s', '--samples',
                    help="excel file with samples in S, F, R, I columns",
                    required=True)

parser.add_argument('-a', '--amplicons',
                    action='store_true',
                    help='Run amplicon workflow instead of shotgun workflow')

parser.add_argument('-q', '--quiet',
                    action='store_false',
                    help='Print out commands')

args = parser.parse_args()

if __name__ == "__main__":
    jak_utils.header()
    files = pd.read_excel(args.samples, index_col=0)

    contam_seqs = jak_utils.get_yaml("contams_db")

    sample_list = []

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        sample_list.append(d)

        print(f'### {d.sample} ###')

        d.verify_read_pairs(echo=args.quiet, run=True)

        if args.amplicons:
            d.contaminant_filtering(contam_seqs, echo=args.quiet, run=True)
        else:
            d.adapter_trimming(contam_seqs, echo=args.quiet, run=True)
            d.contaminant_filtering(contam_seqs, echo=args.quiet, run=True)
            d.quality_filtering(echo=args.quiet, run=True)
