#!/usr/bin/env python

import jak_utils
from jakomics import utilities, colors
from jakomics.utilities import system_call

from jakomics.fastq import FASTQ

import pandas as pd
import argparse
import sys
import os
from tqdm import tqdm

# OPTIONS #####################################################################
parser = argparse.ArgumentParser(
    description='test')

parser.add_argument('-s', '--samples',
                    help="excel file with samples in S, F, R, I columns",
                    required=True)

parser.add_argument('--in_dir',
                    help="Directory with genomes",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual genome files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    default="snippy_out",
                    required=False)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

# FUNCTIONS ###################################################################


def snippy(genome, sample, conda_env="snippy-4.6", snippy_threads=8):

    # create out directory path
    prefix = genome.short_name + "_" + sample.sample
    out_dir = os.path.join(args.out_dir, prefix)

    command = f"snippy --cpus {snippy_threads} --outdir {out_dir} --ref {genome.file_path} --R1 {sample.files[0].file_path} --R2 {sample.files[1].file_path} --prefix {prefix} --cleanup --report"
    #command = f"{command}; conda deactivate"

    snippy_err = system_call(command, echo=False, run=True, return_type='err')

    # write log file to snippy output folder

    out_file = open(os.path.join(out_dir, "snippy_stderr.txt"), 'w')
    for line in snippy_err:
        print(line, file=out_file)
    out_file.close()

# MAIN ########################################################################


if __name__ == "__main__":
    jak_utils.header()
    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    # process fastq files into fastq obbjects
    files = pd.read_excel(args.samples, index_col=0, engine='openpyxl')
    fastq_file_list = []

    for sample, row in files.iterrows():
        d = FASTQ(sample, row)
        fastq_file_list.append(d)

    # get genomes into file objects
    genome_list = utilities.get_files(args.files, args.in_dir, ["gbk", "gbff", "gb", "fa", "fasta"])

    pbar = tqdm(total=len(genome_list) * len(fastq_file_list), desc="Finished", unit=" snippy runs")

    for genome in genome_list:
        for sample in fastq_file_list:
            snippy(genome, sample)
            pbar.update(1)
