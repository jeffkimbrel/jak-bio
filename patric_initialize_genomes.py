import os
import argparse
from tqdm import tqdm
from multiprocessing import Pool
from natsort import natsorted

from jakomics import colors
from jakomics.patric import Patric
import jak_utils

from Bio import SeqIO

# RAST output by default has an incorrect header which leads to a BioPython warning. This will suppress all warnings.
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Run this script on newly downloaded patric genomes', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with top-level patric folders",
                    required=True)

args = parser.parse_args()

# FUNCTIONS ###################################################################

def get_files():
    genome_list = []
    folder_list = [d for d in os.listdir(
        args.in_dir) if os.path.isdir(os.path.join(args.in_dir, d))]
        

    for folder in natsorted(folder_list):
        genome_list.append(Patric(args.in_dir, folder))

    return genome_list

def add_locus_tags(file):

    file.gb_path = os.path.join(args.in_dir, file.full_name, file.full_name + ".gb")
    file.gb_out_path = os.path.join(args.in_dir, file.full_name, file.genome_name + ".gbk")

    if os.path.exists(file.gb_out_path):
        print(f'{colors.bcolors.RED}Overwriting {file.gb_out_path}{colors.bcolors.END}')
        os.remove(file.gb_out_path)

    for seq_record in SeqIO.parse(file.gb_path, "genbank"):
        seq_record.id = seq_record.description

        new_features = []
        for feature in seq_record.features:
            for db_xref in feature.qualifiers['db_xref']:
                if db_xref.startswith('RAST2'):
                    feature.qualifiers['locus_tag'] = db_xref.replace('RAST2:fig|' + file.genome_id, file.genome_name)
                    # print(db_xref, feature.qualifiers['locus_tag'])
            new_features.append(feature)
        seq_record.features = new_features

        output_handle = open(file.gb_out_path, "a")
        SeqIO.write(seq_record, output_handle, "genbank")
        output_handle.close()

# MAIN ########################################################################

if __name__ == "__main__":
    jak_utils.header()

    file_list = get_files()

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(add_locus_tags, file_list), total=len(file_list), desc="Done", unit=" genomes"):
        pass

    pool.close()
