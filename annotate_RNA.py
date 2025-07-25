import argparse
from tqdm import tqdm
from multiprocessing import Pool
import os.path
import sys
import re

import jak_utils
from jakomics import utilities, colors
from jakomics.utilities import system_call

from Bio import SeqIO
from Bio.Seq import Seq

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

# fix out directory
args.out_dir = os.path.abspath(args.out_dir) + '/'
if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)


def getInfo(seq_record):
    
    count = len(str(seq_record.seq))

    return(count)

def run_RFAM(fasta_file):

    genome_size = 0
    for seq_record in SeqIO.parse(fasta_file.file_path, "fasta"):

        count = getInfo(seq_record)
        genome_size += count

    genome_size = genome_size * 2 / 1000000

    rfam_command = f'cmscan -Z {genome_size} --cut_ga --rfam --nohmmonly --tblout {args.out_dir}/{fasta_file.short_name}.tblout --fmt 2 --cpu 1 --clanin {jak_utils.get_yaml("rfam_db")}/Rfam.clanin {jak_utils.get_yaml("rfam_db")}/Rfam.cm {fasta_file.file_path} > {args.out_dir}/{fasta_file.short_name}.cmscan'

    command = 'source activate RFAM && ' + rfam_command + ' && conda deactivate'
    system_call(command, echo=False, run=True)



## MAIN LOOP ###################################################################

if __name__ == "__main__":
    jak_utils.header()

    file_list = utilities.get_files(args.files, args.in_dir, ["fa", "fasta"])

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(run_RFAM, file_list), total=len(file_list), desc="Done", unit=" files"):
        pass

    pool.close()
