import pandas as pd
import argparse
import os
from Bio import SeqIO
import jak_utils
from tqdm import tqdm

jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-m',
                    '--metadata',
                    help="Excel file with BIN and MAG columns for old and new names",
                    required=True)
parser.add_argument('--in_dir',
                    help="Directory with original name bins",
                    required=True)
parser.add_argument('--out_dir',
                    help="Directory to write renamed bins",
                    required=True)

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'
args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Read and Sort Stats File ####################################################
stats = pd.read_excel(args.metadata)
stats['OLD_PATH'] = args.in_dir + stats['BIN']
stats['NEW_PATH'] = args.out_dir + stats['MAG'] + '.fa'

# Interate over bins ##########################################################

pbar = tqdm(total=len(stats.index))

for MAG, row in stats.iterrows():
    counter = 1
    renamed = []
    for record in SeqIO.parse(row['OLD_PATH'], "fasta"):
        record.id = row['MAG'] + '_' + str(counter)
        renamed.append(record)
        counter += 1
    SeqIO.write(renamed, row['NEW_PATH'], "fasta")
    pbar.update(1)
