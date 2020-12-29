from jakomics import colors
import sys
import argparse
import pandas as pd
import jak_utils
jak_utils.header()


# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXXXX')

parser.add_argument('-b', '--bac',
                    help="Path to gtdbtk.bac120.summary.tsv",
                    required=True)

parser.add_argument('-a', '--arc',
                    help="Path to gtdbtk.ar122.summary.tsv",
                    required=True)

parser.add_argument('-o',
                    '--out',
                    help="Path to write output to",
                    required=True)

args = parser.parse_args()

cols = ['user_genome', 'classification', 'fastani_reference', 'fastani_ani',
        'closest_placement_reference', 'closest_placement_ani', 'classification_method', 'note', 'warnings']

df_bac = pd.read_csv(args.bac,
                     sep="\t",
                     index_col=None)

df_bac = df_bac[df_bac.columns.intersection(cols)]

df_arc = pd.read_csv(args.arc,
                     sep="\t",
                     index_col=None)

df_arc = df_arc[df_arc.columns.intersection(cols)]

df = pd.concat([df_bac, df_arc])

df.to_csv(args.out, sep="\t", index=False)
