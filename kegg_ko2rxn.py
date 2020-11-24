import sys
import os
import argparse
from jakomics import colors

import jak_utils
jak_utils.header()

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-k', '--kegg',
                    required=True)

parser.add_argument('-c', '--column',
                    default=1,
                    required=False,
                    type=int,
                    help="Column to replace (0-based)")

args = parser.parse_args()

# kegg db version
v = {}
v_call = 'curl -g -s -S http://rest.kegg.jp/info/kegg'
v_raw = os.popen(v_call).read().split('\n')[1].split("Release ")[1]
print(f'{colors.bcolors.YELLOW}KEGG Release: {v_raw}{colors.bcolors.END}', file=sys.stderr)

## MAKE MAP ###

map = {}
map_call = 'curl -g -s -S http://rest.kegg.jp/link/reaction/orthology/'
map_raw = os.popen(map_call).read()
for line in map_raw.split("\n"):
    split = line.split("\t")

    if len(split) > 1:
        ko = split[0].replace("ko:", "")
        rxn = split[1].replace("rn:", "")

        if ko in map:
            map[ko].append(rxn)
        else:
            map[ko] = [rxn]

## MAP STUFF ##
kegg_raw = [line.strip() for line in open(args.kegg)]
for line in kegg_raw:
    split = line.split("\t")
    if len(split) > 1:

        if split[args.column] in map:
            reactions = map[split[args.column]]
            for rn in reactions:
                print(split[0], rn, sep="\t")
        else:
            print(line)
    else:
        print(line)
