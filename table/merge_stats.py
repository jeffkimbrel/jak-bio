import os
import pandas as pd
from natsort import natsorted

import argparse

# Arguments
parser = argparse.ArgumentParser(description='Merge a folder of stats files, such as from bbtools')

parser.add_argument('-d', '--directory', help="Folder containing the files", required=True)
parser.add_argument('-t', '--text', help="Common file name text to identify files", required=True)
parser.add_argument('-o', '--out', help="Out file name", required=True)

args = parser.parse_args()

all = None

for fileName in natsorted(os.listdir(args.directory)):
    if args.text in fileName:

        split = fileName.split("_")
        file = pd.read_csv(os.path.join(args.directory, fileName), sep="\t")
        file['SAMPLE'] = split[0]

        if all is None:
            all = file

        else:
            all = pd.concat([all, file], sort=False)


all.to_csv(args.out, sep="\t", index=False)
