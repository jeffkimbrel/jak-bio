import sys
import argparse
from collections import Counter # for adding dictionaries

parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-x', '--xref',
    help = "Cross-reference File",
    required = True)

parser.add_argument('-t', '--table',
    help = "Table File",
    required = True)

parser.add_argument('-s', '--split',
    help = "Character to split the xref on, if necessary",
    default = "NA",
    required = False)

parser.add_argument('-c', '--column',
    help = "Column in cross-reference file with new identifier (1-based)",
    default = 2,
    type = int)

parser.add_argument('--divide', '-d',
    action = 'store_true',
    help = 'Divide multiple hits' )

parser.add_argument('--multiples', '-m',
    action = 'store_false',
    help = 'Allow multiples' )

args = parser.parse_args()

args.column -= 1

## CLASSES #####################################################################

xrefs = {}
class XREF:

    '''
    a mapping file between the original and new terms
    '''

    def __init__(self, original):
        self.original = original
        self.new = []

    def addNew(self, new):

        if args.split != "NA":
            newSplit = new.split(args.split)
            for subNew in newSplit:
                self.new.append(subNew)
        else:
            self.new.append(new)

        if args.multiples == True:
            self.new = list(set(self.new))

rows = {}
class ROW:

    '''
    The original table
    '''

    def __init__(self, id, dict):
        self.id = id
        self.dict = dict.copy() # make a copy of the passed dict

    def addAbundance(self, sample, count):
        self.dict[sample] += count

merges = {}
class MERGE:

    def __init__(self, id, dict):
        self.id = id
        self.dict = dict.copy() # make a copy of the passed dict

    def addDictionary(self, rowData, scaleFactor):

        dict = rowData.copy()

        if args.divide:

            for key in dict:
                dict[key] /= scaleFactor

        self.dict = Counter(self.dict) + Counter(dict)

## PROCESS XREF ################################################################

xrefFile = [line.strip() for line in open(args.xref)]

for line in xrefFile:
    split = line.split("\t")
    id = split[0]

    if id not in xrefs:
        xrefs[id] = XREF(id)

    new = "NONE"

    if len(split) > args.column:
        new = split[args.column]

    xrefs[id].addNew(new)

## PROCESS TABLE ###############################################################

tableFile = [line.strip() for line in open(args.table)]

# take the header and convert to an empty dictionary
tableHeader = tableFile.pop(0).split("\t")
tableDict = dict((el, 0) for el in tableHeader)

for line in tableFile:

    split = line.split("\t")

    if (len(tableHeader) + 1) != len(split):
        sys.exit("Uh-oh, does the table header have some extra stuff in the first column?")

    id = split.pop(0)

    # make row class
    rows[id] = ROW(id, tableDict)

    # add abundances
    counter = 0
    while counter < len(tableHeader):
        rows[id].addAbundance(tableHeader[counter], float(split[counter]))
        counter += 1

## MERGE #######################################################################

for original in xrefs:
    for new in xrefs[original].new:
        merges[new] = MERGE(new, tableDict)

for row in rows:
    dict = rows[row].dict
    id = rows[row].id

    if id in xrefs:
        for new in xrefs[id].new:
            merges[new].addDictionary(dict, len(xrefs[id].new))

## PRINT #######################################################################

print(*tableHeader, sep = "\t")

for id in sorted(merges.keys()):

    # get row sum
    sum = 0
    for column in tableHeader:
        sum += merges[id].dict[column]

    if sum > 0:

        print(id, end = "\t")

        for column in tableHeader:
            print(merges[id].dict[column], end = "\t")

        print()
