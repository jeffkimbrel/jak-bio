import sys
import os
import argparse
import uuid
from multiprocessing import Manager, Pool
from tqdm import tqdm

from jakomics import hmm, utilities, colors
from jakomics.table import TABLE, merge_value_counts
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='''
            Find CAZYmes in an amino acid file

            QC_CODE refers to different scoring rules set forth by DBCAN6:
                1A = E-value <= 1e-18 AND HMM coverage >= 35%
                1B = E-value <= 1e-15 AND HMM coverage >= 35%
                2 = E-value <= 1e-5  AND Alignment Length >= 80
                3 = E-value <= 1e-3  AND HMM coverage >= 30%
                0 = Low-quality hit that did not pass any metrics
                ''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--qc',
                    help="QC codes to include",
                    required=False,
                    nargs='*',
                    default=['1A', '1B'])

parser.add_argument('--hmm',
                    help="Write merged HMM results",
                    required=False,
                    default=None)

parser.add_argument('--substrate',
                    help="Write merged substrates results",
                    required=False,
                    default=None)

parser.add_argument('--remove_duplicates',
                    action='store_true',
                    help="""Remove duplicates (hits with the same gene and model)
from the optional merge files. Note, this won't remove
all duplicates (e.g. models such AA1_1 and AA1_2 won't
be considered duplicates). If there is a duplicate, the
row with the highest HMM_COVERAGE will be retained""")

args = parser.parse_args()

manager = Manager()
result_files = manager.list()

# FUNCTIONS ###################################################################


def main(file):

    global result_files

    file.temp_files['temp_log'] = file.id + '.log'
    file.temp_files['temp_out'] = file.id + '.temp.txt'
    file.results_file = file.short_name + '.dbcan8.txt'

    hmm.run_hmmsearch(file.file_path,
                      file.temp_files['temp_log'],
                      file.temp_files['temp_out'],
                      jak_utils.get_yaml("cazyme_db"))

    file.results = hmm.cazymes_to_df(file.temp_files['temp_out'], args.qc)

    if os.path.exists(file.results_file):
        os.remove(file.results_file)
    f = open(file.results_file, 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)

    file.results.to_csv(f, sep="\t", index=False)

    # cleanup
    file.remove_temp()
    result_files.append(file.results_file)


## MAIN LOOP ###################################################################


if __name__ == "__main__":
    jak_utils.header()

    genome_list = utilities.get_files(args.files, args.in_dir, ["faa"])

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(main, genome_list), total=len(genome_list), desc="Finished", unit=" files"):
        pass
    pool.close()

    # # write merged results
    if args.hmm is not None:
        if os.path.exists(args.hmm):
            os.remove(args.hmm)
        f = open(args.hmm, 'a')
        for c in jak_utils.header(r=True):
            print(f'# {c}', file=f)
        print(f'# remove_duplicates={args.remove_duplicates}', file=f)

        print(f'{colors.bcolors.YELLOW}WARNING: {args.hmm} may contain counts of duplicate models per gene{colors.bcolors.END}', file=sys.stderr)
        hmm = merge_value_counts(result_files,
                                 'HMM',
                                 remove_duplicates=args.remove_duplicates)
        hmm.to_csv(f, sep="\t", index=True)

    if args.substrate is not None:
        if os.path.exists(args.substrate):
            os.remove(args.substrate)
        f = open(args.substrate, 'a')
        for c in jak_utils.header(r=True):
            print(f'# {c}', file=f)
        print(f'# remove_duplicates={args.remove_duplicates}', file=f)

        print(f'{colors.bcolors.YELLOW}WARNING: {args.substrate} may contain counts of duplicate models per gene{colors.bcolors.END}', file=sys.stderr)
        substrate = merge_value_counts(result_files,
                                       'SUBSTRATE',
                                       remove_duplicates=args.remove_duplicates)
        substrate.to_csv(f, sep="\t", index=True)
