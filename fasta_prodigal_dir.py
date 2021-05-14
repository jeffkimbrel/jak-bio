import os
from multiprocessing import Pool
import argparse
from tqdm import tqdm

from jakomics import colors, utilities
from jakomics.utilities import system_call

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Run prodigal on a folder of fasta files',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-f', '--files',
                    help="Paths to individual fasta files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--in_dir',
                    help="Directory with nt fasta files",
                    required=True)

parser.add_argument('--out_dir',
                    help="Directory to write .gff and .faa fasta files",
                    required=True)

parser.add_argument('-m',
                    '--meta',
                    help="Enable prodigal meta mode",
                    action='store_true')

args = parser.parse_args()

args.in_dir = os.path.abspath(args.in_dir) + '/'
args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("Creating " + args.out_dir + " because it does not exist\n")
    os.makedirs(args.out_dir)

# FUNCTIONS ###################################################################


def run_prodigal(contig_file):

    command = f"prodigal -q -i {os.path.join(args.in_dir, contig_file.name)} -o {os.path.join(args.out_dir, contig_file.name + '.gff')} -f gff -a {os.path.join(args.out_dir, contig_file.name + '.faa')} -d {os.path.join(args.out_dir, contig_file.name + '.ffn')}"

    if args.meta == True:
        command = command + " -p meta"

    command = 'source activate prodigal-env && ' + command + ' && conda deactivate'
    # systemCall(command, contig_file)
    lines = system_call(command, echo=False, run=True)


if __name__ == "__main__":
    jak_utils.header()

    files = utilities.get_files(args.files, args.in_dir, ['fa', 'fasta'])

    pool = Pool()
    for _ in tqdm(pool.imap_unordered(run_prodigal, files), total=len(files), desc="Finished", unit=" files"):
        pass
    pool.close()
