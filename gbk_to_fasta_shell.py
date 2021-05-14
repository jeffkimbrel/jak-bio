from jakomics import utilities, colors
from jakomics.genome import GENOME

# import gbk_to_fasta
import argparse
import os
import sys

from multiprocessing import Pool
from tqdm import tqdm

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--feature_identifier',
                    help="Feature identifier for genes",
                    required=False,
                    default='locus_tag')

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'


def write_fasta(genome):

    gbk = GENOME(genome)

    # write genes to genomes and gene class dictionary
    gbk.faa_path = os.path.join(args.out_dir, genome.name + ".faa")
    gbk.nt_path = os.path.join(args.out_dir, genome.name + ".ffn")
    gbk.contig_path = os.path.join(args.out_dir, genome.name + ".fa")

    # gbk.genbank_to_fasta(write_faa=gbk.faa_path,
    #     write_nt=gbk.nt_path,
    #     write_contig=gbk.contig_path,
    #     feature_identifier=args.feature_identifier)

    gbk.genbank_to_fasta(write_nt=gbk.nt_path,
                         write_faa=gbk.faa_path,
                         feature_identifier=args.feature_identifier)

# MAIN LOOP ###################################################################


if __name__ == "__main__":
    jak_utils.header()

    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    genome_list = utilities.get_files(args.files, args.in_dir, ["gbk", "gbff", "gb"])

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(write_fasta, genome_list), total=len(genome_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()
