#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
import gzip

from jakomics import colors
from jakomics.fastq import FASTQ, run_info

import jak_utils

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='X')
parser.add_argument('-i1', '--index1', help="Index 1 fastq", required=True)
parser.add_argument('-i2', '--index2', help="Index 2 fastq", required=True)
parser.add_argument('-r1', '--read1', help="Read 1 fastq", required=True)
parser.add_argument('-r2', '--read2', help="Read 2 fastq", required=True)
parser.add_argument('-m',  '--mapping', help="Mapping File", required=True)
parser.add_argument('-b',  '--buffer', help="read buffer size",
                    required=False, type=int, default=2000)
parser.add_argument('--out_dir', help="Output directory", required=False, default="demux")

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)

total = 0

## CLASS #######################################################################

samples = {}


class Sample:
    def __init__(self, sample, barcode):
        self.sample = sample
        self.barcode = barcode
        self.file1 = args.out_dir + sample + "_R1.fastq.gz"
        self.file2 = args.out_dir + sample + "_R2.fastq.gz"
        self.read1 = []
        self.read2 = []
        self.count = 0

        # create empty files
        gzip.open(self.file1, 'w').close()
        gzip.open(self.file2, 'w').close()

    def addPair(self, read1, read2):
        self.read1.append(read1)
        self.read2.append(read2)
        self.count += 1

    def write_to_file(self, final):

        if final == True or len(self.read1) > args.buffer:
            output_handle = gzip.open(self.file1, "at")
            SeqIO.write(self.read1, output_handle, "fastq")
            self.read1 = []
            output_handle.close()

            output_handle = gzip.open(self.file2, "at")
            SeqIO.write(self.read2, output_handle, "fastq")
            self.read2 = []
            output_handle.close()

## FUNCTIONS ###################################################################


def createUndetermined():
    samples['Undetermined'] = Sample('Undetermined', 'Undetermined')


def addSamples():
    mapping_data = pd.read_csv(args.mapping, sep="\t", index_col=0)

    for sample, row in mapping_data.iterrows():
        samples[row['BarcodeSequence']] = Sample(sample, row['BarcodeSequence'])


def readFiles(total_reads):
    total = 0
    errors = 0
    pbar = tqdm(total=total_reads, desc="Processed", unit=" read pairs")

    for read1, read2, index1, index2 in zip(SeqIO.parse(args.read1,  "fastq"),
                                            SeqIO.parse(args.read2,  "fastq"),
                                            SeqIO.parse(args.index1, "fastq"),
                                            SeqIO.parse(args.index2, "fastq")):

        headers_agree = False

        if len(set({read1.id, read2.id, index1.id, index2.id})) == 1:
            headers_agree = True
        else:
            errors += 1

        barcode = index1.seq + index2.seq

        if not barcode in samples:
            barcode = 'Undetermined'

        if headers_agree == True:
            samples[barcode].addPair(read1, read2)

        total += 1
        if total % args.buffer == 0:
            pbar.update(args.buffer)

            if total % (args.buffer * 500) == 0:
                for barcode in samples:
                    samples[barcode].write_to_file(False)

    return(total, errors)


def write_report():
    for sample in samples:
        print(sample.sample, sample.barcode, sample.file1, sample.file2, sample.count, sep="\t")


def main():

    # read count of original
    print(f"{colors.bcolors.GREEN}Getting read count of {args.read1}... {colors.bcolors.END}")
    original = run_info(args.read1)
    for key in original:
        total_reads = original[key]
    print(f"{colors.bcolors.GREEN}...found {total_reads} read pairs {colors.bcolors.END}")

    print(f"{colors.bcolors.GREEN}File buffer size: {args.buffer}, File write cutoff: {args.buffer * 500}{colors.bcolors.END}")

    createUndetermined()
    addSamples()

    total, errors = readFiles(total_reads)

    print(f"Processed {total} total pairs with {errors} pairing errors")

    # final write
    for barcode in samples:
        samples[barcode].write_to_file(True)

    write_report()

## RUN #########################################################################


if __name__ == "__main__":
    jak_utils.header()
    main()
