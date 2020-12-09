#!/usr/bin/env python

import argparse
import pandas as pd

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

parser.add_argument('-m', '--memory',
                    help="Memory to pass to bbtools",
                    required=False,
                    default="Xmx8g")

parser.add_argument('-t', '--threads',
                    help="Threads to pass to bbtools",
                    required=False,
                    default=8)

args = parser.parse_args()


def format_stats(sample, type, stats, sep="|"):
    print(f'{sample}{sep}{type}{sep}IN{sep}{stats["IN"]["READS"]}{sep}{stats["IN"]["BP"]}')
    print(
        f'{sample}{sep}{type}{sep}REMOVED{sep}{stats["REMOVED"]["READS"]}{sep}{stats["REMOVED"]["BP"]}')
    print(f'{sample}{sep}{type}{sep}OUT{sep}{stats["OUT"]["READS"]}{sep}{stats["OUT"]["BP"]}')

    for s in stats:
        if s != "IN":
            if s != "OUT":
                if s != "REMOVED":
                    print(
                        f'{sample}{sep}{type}{sep}{s}{sep}{stats[s]["READS"]}{sep}{stats[s]["BP"]}')


if __name__ == "__main__":
    jak_utils.header()
    files = pd.read_excel(args.samples, index_col=0)

    contam_seqs = jak_utils.get_yaml("contams_db")

    sample_list = []

    print('SAMPLE|STEP|STAT|READS|BP\n---|---|---|---|---')

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        print(d.__dict__)
        sample_list.append(d)

        # print(f'{colors.bcolors.CYAN}\n### {d.sample} ###{colors.bcolors.END}')

        d.verify_read_pairs(echo=args.quiet, run=True)

        if args.amplicons:
            cf = d.contaminant_filtering(contam_seqs, echo=args.quiet,
                                         run=True, mem=args.memory, threads=args.threads)
            format_stats(d.sample, 'CF', cf)

        else:
            rt = d.adapter_trimming(contam_seqs, echo=args.quiet, run=True,
                                    mem=args.memory, threads=args.threads)
            format_stats(d.sample, 'RT', rt)
            cf = d.contaminant_filtering(contam_seqs, echo=args.quiet,
                                         run=True, mem=args.memory, threads=args.threads)
            format_stats(d.sample, 'CF', cf)
            qf = d.quality_filtering(echo=args.quiet, run=True,
                                     mem=args.memory, threads=args.threads)
            format_stats(d.sample, 'QF', qf)
