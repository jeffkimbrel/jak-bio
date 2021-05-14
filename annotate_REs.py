import sys
import os
import argparse
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm
import yaml
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Restriction import RestrictionBatch, Analysis

from Bio.Restriction import AllEnzymes, CommOnly

from jakomics import utilities, colors, file
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

args = parser.parse_args()

# FUNCTIONS ###################################################################


def main(genome):

    original_path = genome.file_path

    if genome.suffix in ['.gb', '.gbk', '.gbff']:
        gbk = GENOME(genome)
        gbk.fa_path = genome.short_name + "_" + genome.id + ".fa"
        gbk.genbank_to_fasta(write_contig=gbk.fa_path)
        genome.file_path = gbk.fa_path
        genome.temp_files['fa'] = gbk.fa_path

    output_file = os.path.join(args.out_dir, genome.name + '.REs.json')

    if os.path.exists(output_file):
        os.remove(output_file)

    for seq_record in SeqIO.parse(genome.file_path, "fasta"):
        a = cleanEnzymeList.search(seq_record.seq)
        # print(a)
        #
        for b in a:
            print(type(b))
            # print(b, len(a[b]), sep="\t")

        #a = Analysis(RestrictionBatch(cleanEnzymeList), seq_record.seq, False)
        # print(type(a))
        with open(output_file, 'w') as outfile:
            json.dump(a, outfile, indent=2)

    #
    # # write to file with comments
    # f = open(output_file, 'a')
    # for c in jak_utils.header(r=True):
    #     print(f'# {c}', file=f)
    # print(f'# FILE: {original_path}', file=f)
    # print(f'# DB: {jak_utils.get_yaml("kofam_db")}', file=f)
    # print(f'# PROFILE: {args.profile}', file=f)
    # print(f'# SCALE: {args.t_scale}', file=f)
    #
    # df.to_csv(f, sep="\t", index=False)

    genome.remove_temp()


# MAIN LOOP ###################################################################

if __name__ == "__main__":
    jak_utils.header()

    # fix out directory
    args.out_dir = os.path.abspath(args.out_dir) + '/'
    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    genome_list = utilities.get_files(args.files, args.in_dir, ["faa", "gbk", "gbff", "gb"])

    if len(genome_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid files were found!{colors.bcolors.END}")

    rawEnzymeList = RestrictionBatch(first=[], suppliers=['N'])
    # rawEnzymeList = RestrictionBatch(first=[], suppliers=['B', 'N', 'K', 'R'])
    print(len(rawEnzymeList))
    cleanEnzymeListNames = []

    for enzyme in rawEnzymeList:
        if enzyme.size > 5:                     # recognition site at least 4 bases
            if enzyme.cut_twice() == False:
                cleanEnzymeListNames.append(enzyme)
    cleanEnzymeList = RestrictionBatch(cleanEnzymeListNames)
    print(len(cleanEnzymeList))

    pool = Pool()
    for _ in tqdm(pool.imap_unordered(main, genome_list), total=len(genome_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()
