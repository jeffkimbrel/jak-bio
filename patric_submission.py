import pandas as pd
import argparse
from tqdm import tqdm

from jakomics.utilities import system_call
from jakomics.file import validate_path
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Submit contigs to Patric. Must log-in to the Patric shell first.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-m', '--metadata',
                    help="Excel with contig information. Required to have CONTIGS_FILE (full path), NAME and TAX_ID columns", required=True)
parser.add_argument('-p', '--patric', help="Patric workspace path", required=True)

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == "__main__":
    jak_utils.header()

    genomes = pd.read_excel(args.metadata, engine='openpyxl')
    pbar = tqdm(total=len(genomes.index), desc="Submitted", unit=" genomes")

    for genome, row in genomes.iterrows():
        validate_path(row['CONTIGS_FILE'])

        command = f"p3-submit-genome-annotation --contigs-file {row['CONTIGS_FILE']} -n \"{row['NAME']}\" -t {row['TAX_ID']} -d Bacteria {args.patric} \"{row['NAME']}\""
        system_call(command, echo=False, run=True)

        pbar.update(1)
