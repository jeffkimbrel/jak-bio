import json
import ast
import sys
import jak_utils
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-c',
                    '--checkm',
                    help="Path to bin_stats_ext.tsv file",
                    required=True)

parser.add_argument('-o',
                    '--out',
                    help="Path to write output to",
                    required=True)

args = parser.parse_args()

#

jak_utils.header()

results = [line.strip() for line in open(args.checkm)]

cols = ["genome", "Completeness", "Contamination", "Genome size", "# scaffolds",
        "Longest scaffold", "N50 (scaffolds)", "marker lineage", "GC", "Coding density", "# markers"]

df = pd.DataFrame(columns=cols)

for line in results:
    sample, js = line.split("\t")
    # sample = sample + ".fa"
    js = ast.literal_eval(js)
    js['genome'] = sample

    series = pd.Series(js)
    df = df.append(series, ignore_index=True)

df = df[df.columns.intersection(cols)]

df.to_csv(args.out, sep="\t", index=False)
