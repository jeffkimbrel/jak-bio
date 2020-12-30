import sys
import os
import argparse
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm
import yaml

from jakomics import utilities, colors, kegg, file
from jakomics.genome import GENOME
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
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

parser.add_argument('--profile',
                    help="kofamscan profile",
                    default='prokaryote',
                    required=False)

parser.add_argument('--t_scale',
                    help="Threshold Scale",
                    default=1.0,
                    type=float,
                    required=False)

args = parser.parse_args()

# FUNCTIONS ###################################################################


def main(genome):

    original_path = genome.file_path

    if genome.suffix in ['.gb', '.gbk', '.gbff']:
        gbk = GENOME(genome)
        gbk.faa_path = genome.short_name + "_" + genome.id + ".faa"
        gbk.genbank_to_fasta(write_faa=gbk.faa_path)
        genome.file_path = gbk.faa_path
        genome.temp_files['faa'] = gbk.faa_path

    output_file = os.path.join(args.out_dir, genome.name + '.KofamKOALA.txt')

    if os.path.exists(output_file):
        os.remove(output_file)

    hits = kegg.run_kofam(genome.file_path,
                          args.profile,
                          os.path.join(jak_utils.get_yaml("kofam_db"), 'ko_list'),
                          t_scale=args.t_scale)

    df = kegg.kofam_to_df(hits)

    # write to file with comments
    f = open(output_file, 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)
    print(f'# FILE: {original_path}', file=f)
    print(f'# DB: {jak_utils.get_yaml("kofam_db")}', file=f)
    print(f'# PROFILE: {args.profile}', file=f)
    print(f'# SCALE: {args.t_scale}', file=f)

    df.to_csv(f, sep="\t", index=False)

    genome.remove_temp()


# MAIN LOOP ###################################################################

if __name__ == "__main__":
    jak_utils.header()

    # fix out directory
    args.out_dir = os.path.abspath(args.out_dir) + '/'
    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    # fix and validate hal file path
    if args.profile == 'prokaryote':
        args.profile = os.path.join(jak_utils.get_yaml("kofam_db"), 'prokaryote.hal')
    file.validate_path(args.profile)

    genome_list = utilities.get_files(args.files, args.in_dir, ["faa", "gbk", "gbff", "gb"])

    if len(genome_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid files were found!{colors.bcolors.END}")

    pool = Pool()
    for _ in tqdm(pool.imap_unordered(main, genome_list), total=len(genome_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()
