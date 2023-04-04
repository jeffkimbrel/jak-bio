#!/usr/bin/env python

import argparse
import os
from tqdm import tqdm
import multiprocessing

from Bio import SeqIO
import gzip

import jak_utils
from jakomics import colors, utilities

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Takes a folder of fastq files and splits each file into subfiles based on the run IDs of the fastq headers.')

parser.add_argument('--in_dir',
                    help="Directory with fastq files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual fastq files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Parent to write child directories and fastq files to",
                    default="fastq_split",
                    required=False)

parser.add_argument('--lane_split',
                    action='store_true',
                    help='Samples run on multiple lanes will be split')

args = parser.parse_args()

# FUNCTIONS

def split_file(file):
    file.runID = {}
    with gzip.open(file.file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            split = record.id.split(":")

            # make the run ID from the header. The delimiter is normally a ":", but "_" is used here so it is directory name compatible
            if args.lane_split:
                runID = f"{split[0]}_{split[1]}_{split[2]}_{split[3]}"
            else:
                runID = f"{split[0]}_{split[1]}_{split[2]}"
    
            if runID in file.runID:
                file.runID[runID].append(record)
            else:
                file.runID[runID] = [record]
    
    for runID in file.runID:
        runID_path = os.path.join(args.out_dir, runID)
        
        if not os.path.exists(runID_path):
            os.makedirs(runID_path)

        output_handle = gzip.open(os.path.join(runID_path, file.name), "at")
        SeqIO.write(file.runID[runID], output_handle, "fastq")
        output_handle.close()



#

if __name__ == "__main__":
    fastq_files = utilities.get_files(args.files, args.in_dir, ['fq.gz', 'fastq.gz'])
    
    jak_utils.header()

    # make out directory
    args.out_dir = os.path.abspath(args.out_dir) + '/'
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # run
    pool = multiprocessing.Pool()
    for _ in tqdm(pool.imap_unordered(split_file, fastq_files), total=len(fastq_files), desc="Finished", unit=" fastq files"):
        pass
    pool.close()
