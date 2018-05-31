# given a taxa table, this will added "unclassified_XXX" to NA

import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-t', '--tax',
    help = "Tax Table",
    required = True)

args = parser.parse_args()

###

with open(args.tax) as f:
    for line in f:
        line = line.strip()
        split = line.split("\t")

        print(split[0], end = "\t")

        counter = 1
        lastGood = ""
        while counter < len(split):

            # if good
            if split[counter] != "NA":
                lastGood = split[counter]
                print(split[counter], end = "\t")
                counter += 1
            else:
                print("unclassified_" + lastGood, end = "\t")
                counter += 1

        print()
