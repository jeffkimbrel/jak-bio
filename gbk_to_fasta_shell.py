from jakomics import utilities, colors
import gbk_to_fasta
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

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'


def write_fasta(genome):

    # write genes to genomes and gene class dictionary
    genome.faa_path = os.path.join(args.out_dir, genome.name + ".faa")
    genome.nt_path = os.path.join(args.out_dir, genome.name + ".ffn")
    genome.contig_path = os.path.join(args.out_dir, genome.name + ".fa")

    genome.genes = gbk_to_fasta.main(genome.file_path,
                                     write_faa=genome.faa_path,
                                     write_nt=genome.nt_path,
                                     write_contig=genome.contig_path,
                                     return_gene_dict=False)


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
