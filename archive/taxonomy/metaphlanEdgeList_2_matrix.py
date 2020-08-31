import os, sys

import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-f', '--folder',
    help = "Folder path with edgelist files",
    required = True)

parser.add_argument('-l', '--level',
    help = "Taxonomic Level (kpcofgst)",
    default = "s")

args = parser.parse_args()

args.level = args.level.lower()+"__"

# dictionary

metaDict = {}
fileList = []
taxaList = []

# Process dem files

path = args.folder
dirs = os.listdir( path )

for fileName in dirs:
    if fileName.endswith('edgeList.txt'):
        filePath = path + "/" + fileName
        data = open(filePath, 'rt')
        (name,ext) = os.path.splitext(fileName)
        name = name.split(".")[0]

        if name not in fileList:
            fileList.append(name)

        if name not in metaDict:
            metaDict[name] = {}

        for line in data:
            line = line.strip()
            A,B,C = line.split("\t")

            if B.startswith(args.level):

                if B not in taxaList:
                    taxaList.append(B)

                if B in metaDict[name]:
                    metaDict[name][B] += C
                else:
                    metaDict[name][B] = C

        data.close()

# Print header

for sample in fileList:
    print(sample, end="\t")
print()

for taxa in taxaList:
    print(taxa, end = "\t")

    for sample in fileList:
        if taxa in metaDict[sample]:
            print(metaDict[sample][taxa], end = "\t")
        else:
            print(0, end = "\t")

    print()
