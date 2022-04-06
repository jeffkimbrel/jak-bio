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

parser.add_argument('--card_db',
                    help="card.json file",
                    required=False,
                    default="/Users/kimbrel1/Dropbox/Lab/Resources/CARD/card.json")

parser.add_argument('--include_loose',
                    action='store_true',
                    help="Include loose hits")

args = parser.parse_args()

def rgi(file):
    # make temp file without asterisks
    file.temp_files['asterisk_free'] = f'{file.file_path}.temp'

    seqs = []
    for seq_record in SeqIO.parse(file.file_path, "fasta"):
        if seq_record.seq.endswith("*"):
            seq_record.seq = Seq(re.sub(r"\*$", "", str(seq_record.seq)))
        seqs.append(seq_record)

    SeqIO.write(seqs, file.temp_files['asterisk_free'], "fasta")

    command = f"rgi main --input_sequence {file.temp_files['asterisk_free']} --output_file {file.short_name}.CARD --input_type protein --alignment_tool diamond --local --clean -n 1"
    if args.include_loose:
        command = command + " --include_loose"
    
    rgi_err = system_call(command, echo=False, run=True, return_type='err')
    
    file.remove_temp()

## MAIN LOOP ###################################################################

if __name__ == "__main__":
    jak_utils.header()

    file_list = utilities.get_files(args.files, args.in_dir, ["faa"])

    # load CARD database
    if os.path.isfile(args.card_db):
        if os.path.isdir('localDB'):
            print(f"{colors.bcolors.YELLOW}CARD localDB already exists... using that copy{colors.bcolors.END}")
        else:
            print(f"{colors.bcolors.YELLOW}Loading card.json from {args.card_db}{colors.bcolors.END}")
            system_call('rgi load --card_json ~/Dropbox/Lab/Resources/CARD/card.json --local', echo=False, run=True, return_type='out')
            print(f"{colors.bcolors.YELLOW}WARNING: Sometimes the localDB fails the first time it is run... if you get errors, just cancel and start again.{colors.bcolors.END}")
        
        print(f"{colors.bcolors.YELLOW}CARD DB Version: {system_call('rgi database --version', echo=False, run=True, return_type='out')[0]}{colors.bcolors.END}")
        print(f"{colors.bcolors.YELLOW}RGI Version: {system_call('rgi main --version', echo=False, run=True, return_type='out')[0]}{colors.bcolors.END}")
        
    else:
        print(f"{colors.bcolors.RED}NO FILE FOUND AT {args.card_db}{colors.bcolors.END}")
        sys.exit()


    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(rgi, file_list), total=len(file_list), desc="Done", unit=" files"):
        pass

    pool.close()
