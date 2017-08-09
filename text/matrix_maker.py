import sys
import os
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-d', '--directory', help = "Folder containing the files", required = True)
parser.add_argument('-n', '--name', help = "Common feature of name to identify files", required = True)
parser.add_argument('-c', '--column', default = 3, help = "Column with counts to be used" )
parser.add_argument('-m', '--minimum', default = 0, help = "Minimum value" )

args = parser.parse_args()
args.column = int(args.column) - 1
args.minimum = int(args.minimum)


# global variables
dirs = os.listdir(args.directory)
dictionary = {}
columnNames = []
#rowNames = []

for fileName in dirs:
    if args.name in fileName:

        # Get Sample Name
        fileNameSplit = fileName.split(".")
        sampleName = fileNameSplit[0]
        columnNames.append(sampleName)

        # Process file
        with open(args.directory+"/"+fileName) as f:
            for line in f:
                split = line.strip().split("\t")

                rowName = split[0]
                #if rowName not in rowNames:
                #    rowNames.append(rowName)

                value = int(split[args.column])
                if value > args.minimum:
                    if rowName in dictionary:
                        dictionary[rowName][sampleName] = value
                    else:
                        dictionary[rowName] = {sampleName : value}


# print header
for sampleName in columnNames:
    print(sampleName, end = "\t")
print()

for rowName in sorted(dictionary.keys()):
    print(rowName, end = "\t")
    for sampleName in columnNames:
        if sampleName in dictionary[rowName]:
            print(dictionary[rowName][sampleName], end = "\t")
        else:
            print(0, end = "\t")
    print()
