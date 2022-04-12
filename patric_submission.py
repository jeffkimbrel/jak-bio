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
                    help="Excel/tsv with contig information. Required to have CONTIGS_FILE (full path), GENOME and NCBI_ID columns", required=True)
parser.add_argument('-p', '--patric', help = "Patric workspace path", required = True)

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == "__main__":
    jak_utils.header()

    try:
        genomes = pd.read_excel(args.metadata, engine='openpyxl')
    except:
        try:
           genomes = pd.read_csv(args.metadata, sep = "\t")
           
        except:
            print(f"Give me some usable file.")




    
    pbar = tqdm(total=len(genomes.index), desc="Submitted", unit=" genomes")

    for genome, row in genomes.iterrows():
        validate_path(row['CONTIGS_FILE'])

        command = f"p3-submit-genome-annotation --contigs-file \"{row['CONTIGS_FILE']}\" -n \"{row['NAME']}\" -t {row['NCBI_ID']} -d {row['TYPE']} {args.patric} \"{row['GENOME']}\""
        system_call(command, echo=False, run=True)

        pbar.update(1)
