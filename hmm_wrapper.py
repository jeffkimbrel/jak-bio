import os
import argparse
from multiprocessing import Manager, Pool
from tqdm import tqdm

from jakomics import hmm, utilities, colors, file
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

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

parser.add_argument('-d', '--database',
                    help="HMM database",
                    required=True)

args = parser.parse_args()

# manager = Manager()
# result_files = manager.list()

hmm_file = file.FILE(args.database)

# FUNCTIONS ###################################################################


def main(file):

    global hmm_file

    file.temp_files['temp_log'] = file.id + '.log'
    file.results_file = file.short_name + '.hmm.txt'

    hmm.run_hmmsearch(file.file_path,
                      file.temp_files['temp_log'],
                      os.path.join(args.out_dir, file.short_name + '.' + hmm_file.short_name + '.txt'),
                      hmm_file.file_path,
                      echo=False, run=True)

    # cleanup
    file.remove_temp()


## MAIN LOOP ###################################################################

if __name__ == "__main__":
    jak_utils.header()

    # make output directory if it doesn't exist
    args.out_dir = os.path.abspath(args.out_dir) + '/'
    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)    
    
    genome_list = utilities.get_files(args.files, args.in_dir, ["faa", "feature_protein.fasta"])

    pool = Pool()
    for _ in tqdm(pool.imap_unordered(main, genome_list), total=len(genome_list), desc="Finished", unit=" files"):
        pass
    pool.close()