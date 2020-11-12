import os
from multiprocessing import Pool
import argparse
import subprocess
from natsort import natsorted
from tqdm import tqdm

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Run prodigal on a folder of fasta files',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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


def get_files():
    files = []
    dirs = os.listdir(args.in_dir)
    for fileName in dirs:
        if fileName.endswith('fa'):
            files.append(fileName)

    return natsorted(files)


def systemCall(command, contig_file):
    # print("\n>$ "+command)
    subprocess.run('source activate ~/bin/prodigal/ && ' +
                   command + ' && conda deactivate', shell=True)

    # print("\nDone: " + contig_file)


def run_prodigal(contig_file):

    command = 'prodigal -q -i ' + args.in_dir + contig_file + ' -o ' + args.out_dir + contig_file + \
        '.gff -f gff -a ' + args.out_dir + contig_file + '.faa -d ' + args.out_dir + contig_file + '.ffn'

    if args.meta == True:
        command = command + " -p meta"

    systemCall(command, contig_file)


if __name__ == "__main__":

    files = get_files()
    #
    # pool = Pool()
    # pool.map(run_prodigal, files)
    # pool.close()
    # pool.join()

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(run_prodigal, files), total=len(files), desc="Finished", unit=" files"):
        pass
    pool.close()
