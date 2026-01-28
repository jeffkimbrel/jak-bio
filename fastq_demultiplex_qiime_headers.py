#!/usr/bin/env python

import argparse
import os
import sys
from Bio import SeqIO
from tqdm import tqdm
import gzip

from jakomics import colors
from jakomics.fastq import FASTQ, run_info

import jak_utils

# import warnings
# from Bio import BiopythonWarning
# warnings.simplefilter('ignore', BiopythonWarning)


## OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Demultiplex de-interleaved fastq files with de-interleaved barcode files and a sample sheet')

parser.add_argument('-r1', '--read1', help="Read 1 fastq", required=True)
parser.add_argument('-r2', '--read2', help="Read 2 fastq", required=True)

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
    def __init__(self, sample):
        self.sample = sample
        self.file1 = args.out_dir + sample + "_R1.fastq.gz"
        self.file2 = args.out_dir + sample + "_R2.fastq.gz"
        self.read1 = []
        self.read2 = []
        self.singles = {}
        self.count = 0
        self.errors = 0
        self.fixed = 0

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

def readFiles(total_reads):
    total = 0
    errors = 0
    pbar = tqdm(total=total_reads, desc="Processed", unit=" read pairs")

    for read1, read2 in zip(SeqIO.parse(args.read1,  "fastq"),
                            SeqIO.parse(args.read2,  "fastq")):
    
        sample1, pair1, direction1 = parse_header(read1)
        sample2, pair2, direction2 = parse_header(read2)

        

        if not sample1 in samples:
            samples[sample1] = Sample(sample1)

        if not sample2 in samples:
            samples[sample2] = Sample(sample2)
   
        if len(set({pair1, pair2})) == 1:
            # overwrite existing header with sample_pair part removed
            read1.description = " ".join(read1.description.split(" ")[1:])
            read2.description = " ".join(read2.description.split(" ")[1:])
            
            samples[sample1].addPair(read1, read2)
        else:
            errors += 1
            samples[sample1].errors += 1

            addSingles(read1, read2)
        
        total += 1

        if total % args.buffer == 0:
            pbar.update(args.buffer)

            if total % (args.buffer * 500) == 0:
                for sample in samples:
                    samples[sample].write_to_file(False)
                    

    return(total, errors)

def addSingles(read1, read2):
    sample1, pair1, direction1 = parse_header(read1)
    sample2, pair2, direction2 = parse_header(read2)

    # overwrite existing header with sample_pair part removed
    read1.description = " ".join(read1.description.split(" ")[1:])
    read2.description = " ".join(read2.description.split(" ")[1:])
    
    if pair1 in samples[sample1].singles:
        samples[sample1].singles[pair1][direction1] = read1
    else:
        samples[sample1].singles[pair1] = {}
        samples[sample1].singles[pair1][direction1] = read1

    if pair2 in samples[sample2].singles:
        samples[sample2].singles[pair2][direction2] = read2
    else:
        samples[sample2].singles[pair2] = {}
        samples[sample2].singles[pair2][direction2] = read2
    
def parse_header(read):

    # returns the sample, pair id and direction (0 or 1)

    header = read.description.split(" ")
    return(header[0].split("_")[0], header[1], int(header[2].split(":")[0]) - 1)

def write_report():
    for sample in samples:
        print(samples[sample].sample, samples[sample].file1, samples[sample].file2, samples[sample].count, samples[sample].errors, samples[sample].fixed, sep="\t")

def get_read_count(file, print_results = True):

    original = run_info(file)
    total_reads = 0
    for key in original:
        total_reads += original[key]
    
    if print_results:
        print(f"{colors.bcolors.GREEN}Found {total_reads} reads in {file} {colors.bcolors.END}")

    return(total_reads)

def main():

    # read count of original
    
    total_reads = get_read_count(args.read1, print_results = True)

    print(f"{colors.bcolors.GREEN}File buffer size: {args.buffer}, File write cutoff: {args.buffer * 500}{colors.bcolors.END}")
    
    total, errors = readFiles(total_reads)
    fixed = 0
    # for sample in samples:
    #     print(sample, samples[sample].count, samples[sample].read1)

    print(f"Processed {total} total pairs with {errors} pairing errors")

    # final write
    for sample in samples:
        # print(sample, len(samples[sample].singles))
        for single in samples[sample].singles:
            if len(samples[sample].singles[single]) == 1:
                print(samples[sample].singles[single])
            else:
                samples[sample].addPair(samples[sample].singles[single][0], samples[sample].singles[single][1])
                samples[sample].fixed += 1
                fixed += 1
            # for direction in samples[sample].singles[single]:
            #     print(direction, samples[sample].singles[single][direction].description)
        # print(samples[sample].singles)

        samples[sample].write_to_file(True)

    print(f"Fixed {fixed} total pairs with {errors}")

    write_report()

## RUN #########################################################################


if __name__ == "__main__":
    jak_utils.header()
    main()
