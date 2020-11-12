from jakomics import colors
import sys
import argparse
import pandas as pd
import jak_utils
jak_utils.header()

print(
    f'{colors.bcolors.YELLOW}NOTE: this script is currently not doing anything with the gtdbtk.ar122.summary.tsv file if given.{colors.bcolors.END}', file=sys.stderr)

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXXXX')

parser.add_argument('-b', '--bac',
                    help="Path to gtdbtk.bac120.summary.tsv",
                    required=True)

parser.add_argument('-a', '--ar',
                    help="Path to gtdbtk.ar122.summary.tsv",
                    required=False)

parser.add_argument('-o',
                    '--out',
                    help="Path to write output to",
                    required=True)

args = parser.parse_args()

cols = ['user_genome', 'classification', 'fastani_reference', 'fastani_ani',
        'closest_placement_reference', 'closest_placement_ani', 'classification_method', 'note', 'warnings']

df = pd.read_csv(args.bac,
                 sep="\t",
                 index_col=None)

df = df[df.columns.intersection(cols)]
df.to_csv(args.out, sep="\t", index=False)
