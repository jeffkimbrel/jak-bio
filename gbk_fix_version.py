from jakomics import utilities, colors
from jakomics.genome import GENOME, increment_gbk_version, add_gbk_comment

# import gbk_to_fasta
import argparse
import os

from multiprocessing import Pool
from tqdm import tqdm
from Bio import SeqIO

import jak_utils

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

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

args = parser.parse_args()
# args.out_dir = os.path.abspath(args.out_dir) + '/'

def fix_version(gbk):

    gbk.out = f"{gbk.dir}/{gbk.short_name}.1.gbk"

    for seq_record in SeqIO.parse(gbk.file_path, "genbank"):

        seq_record.annotations['molecule_type'] = seq_record.annotations['molecule_type'].upper()
        seq_record = increment_gbk_version(seq_record)
        seq_record = add_gbk_comment(seq_record, "version increment to make KBase compatible")

        output_handle = open(gbk.out, "a")
        SeqIO.write(seq_record, output_handle, "genbank")
        output_handle.close()

# MAIN LOOP ###################################################################


if __name__ == "__main__":
    jak_utils.header()
    # if not os.path.exists(args.out_dir):
    #     print("\nCreating directory " + args.out_dir)
    #     os.makedirs(args.out_dir)

    genome_list = utilities.get_files(args.files, args.in_dir, ["gbk", "gbff", "gb"])

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(fix_version, genome_list), total=len(genome_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()
