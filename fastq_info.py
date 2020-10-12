#!/usr/bin/env python

import argparse
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

from jakomics.fastq import FASTQ, run_info
from jakomics import colors

import jak_utils


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


def get_info(sample):
    for fastq_file in sample.files:
        if args.md5:
            fastq_file.get_md5()
        else:
            fastq_file.md5 = "NA"

        run_info_results = run_info(fastq_file.file_path)

        if len(run_info_results) > 1:
            print(
                f"\n{colors.bcolors.RED}WARNING: {fastq_file.name} appears to be a combination of different Illumina runs!{colors.bcolors.END}")

        for result in run_info_results:
            print(sample.sample, fastq_file.name, sample.type, fastq_file.read,
                  result, run_info_results[result], fastq_file.md5, sep="\t")


if __name__ == "__main__":
    jak_utils.header()
    files = pd.read_excel(args.samples, index_col=0)

    sample_list = []

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        sample_list.append(d)

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(get_info, sample_list), total=len(sample_list), desc="Finished", unit=" samples"):
        pass
    pool.close()
