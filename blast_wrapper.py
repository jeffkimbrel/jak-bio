import argparse
import os

from multiprocessing import Pool
from tqdm import tqdm
import pandas as pd

from jakomics import utilities, blast, colors
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
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

parser.add_argument('--blast_db', '-db',
                    help="Path to blast database/file",
                    required=True)

parser.add_argument('--blast_type', '-t',
                    help="Select either blastn or blastp",
                    required=True)

args = parser.parse_args()

## 

if args.blast_type == 'blastn':
    blast_type = "nucl"
elif args.blast_type == "blastp":
    blast_type = "prot"
else:
    sys.exit("Select either blastn or blastp with the -t flag")


def blast_files(file):
    blast_result = blast.run_blast(type=blast_type,
                q=file.file_path,
                db=args.blast_db,
                e=1e-7,
                make=False)

    blast_results = blast.blast_to_df(blast_result)
    file.results_file = os.path.join(args.out_dir, file.name + '.blastp.txt')

    if os.path.exists(file.results_file):
        os.remove(file.results_file)
    f = open(file.results_file, 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)

    blast_results.to_csv(f, sep="\t", index=False)

if __name__ == "__main__":
    jak_utils.header()

    args.out_dir = os.path.abspath(args.out_dir) + '/'
    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    file_list = utilities.get_files(args.files, args.in_dir, ["fasta", "ffn", "faa", "fa"])

    blast.make_blast_db(blast_type, args.blast_db)

    if len(file_list) > 8:
        pool = Pool(processes=8)
    else:
        pool = Pool(processes=len(file_list))

    for _ in tqdm(pool.imap_unordered(blast_files, file_list), total=len(file_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()