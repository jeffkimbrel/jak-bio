from natsort import natsorted
import argparse
import os
from tqdm import tqdm
from multiprocessing import Pool

import jak_utils
from jakomics.utilities import system_call
from jakomics import colors

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Download genomes in a Patrick workspace. Must log-in to the Patric shell first.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-p', '--patric', help="Patric workspace path", required=True)
parser.add_argument('--out_dir', help="Directory to download the Patric data to", required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

# Functions ###################################################################


def get_genome_list():
    genome_list = []

    command = f'p3-ls -l --type {args.patric}'
    # p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = p.communicate()
    ls_results = system_call(command, echo=False, return_type="out")

    for result in ls_results:
        split = result.split()
        if len(split) > 0:
            if split[6] == 'job_result':
                genome_list.append('\ '.join(split[7:]))

    if len(genome_list) == 0:
        raise Exception(
            f"{colors.bcolors.RED}No genomes found... check path (-p){colors.bcolors.END}")
    else:
        print(f"{colors.bcolors.GREEN}Found {len(genome_list)} genome directories at {args.patric}{colors.bcolors.END}")
        return genome_list


def download_patric_folder(genome):

    command_main = f'p3-cp -R ws:{args.patric}/.{genome} {os.path.join(args.out_dir, genome)}'

    system_call(command_main, echo=False, run=True)

# MAIN ########################################################################


if __name__ == "__main__":
    jak_utils.header()

    genomes = get_genome_list()
  
    pool = Pool(processes=4)
    for _ in tqdm(pool.imap_unordered(download_patric_folder, genomes), total=len(genomes), desc="Downloaded", unit=" genomes"):
        pass

    pool.close()
