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

parser.add_argument('-w', '--workflow',
                    help="amplicon (a) or shotgun (s) workflow",
                    required=False,
                    default="s")


args = parser.parse_args()

if __name__ == "__main__":
    jak_utils.header()
    files = pd.read_excel(args.samples, index_col=0)

    contam_seqs = jak_utils.get_yaml("contams_db")

    sample_list = []

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        sample_list.append(d)

        d.current_name = d.sample
        d.verify_read_pairs()

        if args.workflow == 's':
            d.adapter_trimming(contam_seqs, echo=True, run=True)

        if args.workflow == 'a':
            d.contaminant_filtering(contam_seqs, echo=True, run=True)
