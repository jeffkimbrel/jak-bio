#!/usr/bin/env python

import jak_utils
from jakomics import colors
from jakomics.fastq import FASTQ, run_info
import argparse
import os
import pandas as pd
from tqdm import tqdm
import multiprocessing
multiprocessing.set_start_method("fork")  # python 3.8 fix


# OPTIONS #####################################################################
parser = argparse.ArgumentParser(
    description='test')

parser.add_argument('-s', '--samples',
                    help="excel file with samples in S, F, R, I columns",
                    required=True)

parser.add_argument('--md5',
                    action='store_true',
                    help='Run md5 (very slow)')

parser.add_argument('--out',
                    help="File to write results to",
                    default="fastq_info_results.txt",
                    required=False)

args = parser.parse_args()


def get_info(sample):

    for fastq_file in sample.files:
        if args.md5:
            fastq_file.get_md5()
        else:
            fastq_file.md5 = "NA"

        run_info_results = run_info(fastq_file.file_path)

        file_results = {"SAMPLE": sample.sample,
                        "FILE": fastq_file.file_path,
                        "MD5": fastq_file.md5,
                        "TYPE": sample.type,
                        "PAIR": fastq_file.read,
                        "TOTAL_READS": 0,
                        "RUN_INFO": {},
                        }

        for result in run_info_results:
            file_results["RUN_INFO"][result] = run_info_results[result]
            file_results["TOTAL_READS"] += run_info_results[result]

        results[f'{sample.sample}_{fastq_file.read}'] = file_results


if __name__ == "__main__":
    manager = multiprocessing.Manager()
    results = manager.dict()

    jak_utils.header()
    files = pd.read_excel(args.samples, index_col=0, engine='openpyxl')

    sample_list = []

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        sample_list.append(d)

    pool = multiprocessing.Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(get_info, sample_list), total=len(sample_list), desc="Finished", unit=" samples"):
        pass
    pool.close()

    df = pd.DataFrame(columns=['INDEX', 'SAMPLE', 'FILE', 'MD5',
                               'TYPE', 'PAIR', 'TOTAL_READS', 'RUN_INFO'])

    for result in results:
        df = df.append(
            pd.Series(data={
                'INDEX': result,
                'SAMPLE': results[result]["SAMPLE"],
                'FILE': results[result]["FILE"],
                'MD5': results[result]["MD5"],
                'TYPE': results[result]["TYPE"],
                'PAIR': results[result]["PAIR"],
                'TOTAL_READS': results[result]["TOTAL_READS"],
                'RUN_INFO': str(results[result]["RUN_INFO"]),

            }
            ),
            ignore_index=True)

    df = df.sort_values(by=['INDEX'])

    # write to file with comments
    if os.path.exists(args.out):
        os.remove(args.out)
    f = open(args.out, 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)
    for arg in vars(args):
        print(f'# ARG {arg} = {getattr(args, arg)}', file=f)

    df.to_csv(f, sep="\t", index=False)
