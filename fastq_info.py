#!/usr/bin/env python

import argparse
import pandas as pd

from jakomics.fastq import FASTQ, run_info
from jakomics import colors

import jak_utils
jak_utils.header()

# OPTIONS #####################################################################
parser = argparse.ArgumentParser(
    description='test')

parser.add_argument('-s', '--samples',
                    help="excel file with samples in S, F, R, I columns",
                    required=False)

parser.add_argument('--md5',
                    action='store_true',
                    help='Run md5 (very slow)')

args = parser.parse_args()

files = pd.read_excel(args.samples, index_col=0)

for sample, row in files.iterrows():
    print(f"Working on {sample}")
    d = FASTQ(sample, row)
    for fastq_file in d.files:
        if args.md5:
            fastq_file.get_md5()

        run_info_results = run_info(fastq_file.file_path)

        if len(run_info_results) > 1:
            print(
                f"\n{colors.bcolors.RED}WARNING: {fastq_file.name} appears to be a combination of different Illumina runs!{colors.bcolors.END}")

        for result in run_info_results:
            print(d.sample, fastq_file.name, d.type, fastq_file.read,
                  result, run_info_results[result], sep="\t")
