import argparse

import pandas as pd
import os
import sys

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-m',
                    '--metadata',
                    help="Metadata file",
                    required=True)
parser.add_argument(
    '--out_dir',
    help="Directory to move the files to",
    required=True)
parser.add_argument('--filetype',
                    help="File type",
                    default="faa",
                    required=False)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'


if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Read and Sort Stats File ####################################################
genomes = pd.read_csv(args.metadata, sep="\t")

for genome, row in genomes.iterrows():
    name = str(row['organism_name'])
    if pd.notna(row['infraspecific_name']):
        name = name + " " + str(row['infraspecific_name'])
    name = name.replace(" ", "_")
    name = name.replace("strain=", "")

    name = name.replace("(", "")
    name = name.replace(")", "")
    name = name.replace("*", "_")
    name = name.replace("/", "_")

    genomes.loc[genome, 'NAME'] = name

# check for duplicates
if genomes['NAME'].duplicated().any():
    print(genomes['NAME'].duplicated() == True)
    sys.exit("Ooops, looks like there are some duplicate genomes")

genomes = genomes.set_index('NAME')

for name, row in genomes.iterrows():
    print(name, row['local_filename'], sep="\t")
    command = "gunzip -qkc " + row['local_filename'] + \
        " > " + args.out_dir + name + '.' + args.filetype
    os.system(command)
