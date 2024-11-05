import pandas as pd
import argparse
import os
from Bio import SeqIO
from tqdm import tqdm

from jakomics import colors, file
import jak_utils
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
parser.add_argument('--extension',
                    help="File extension",
                    default = "fa",
                    required=False)

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'
args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print(f"{colors.bcolors.BLUE}Creating directory {args.out_dir}{colors.bcolors.END}")
    os.makedirs(args.out_dir)

# Read and Sort Stats File ####################################################
stats = pd.read_excel(args.metadata, engine='openpyxl')
stats['NEW_PATH'] = f"{args.out_dir}{stats['MAG']}.{args.extension}"

# Interate over bins ##########################################################

pbar = tqdm(total=len(stats.index))

for MAG, row in stats.iterrows():
    bin_path = f"{args.in_dir}{row['BIN']}.{args.extension}"
    mag_path = f"{args.out_dir}{row['MAG']}.fa"

    f = file.FILE(bin_path)
    if f.check_files_exist(exit_if_false = False):
        
        counter = 1
        renamed = []
        for record in SeqIO.parse(f.file_path, "fasta"):
            record.id = row['MAG'] + '_' + str(counter)
            renamed.append(record)
            counter += 1
        SeqIO.write(renamed, mag_path, "fasta")
    pbar.update(1)