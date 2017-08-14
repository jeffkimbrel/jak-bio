import sys
import os
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-d', '--directory', help = "Folder containing the files", required = True)
parser.add_argument('-n', '--name', help = "Common feature of name to identify files", required = True)
parser.add_argument('-c', '--column', default = 3, type = int, help = "Column with counts to be used (1-based)" )
parser.add_argument('-m', '--minimum', default = 0, type = int, help = "Minimum value" )
parser.add_argument('-o', '--out', help = "Out file name", required = True)

args = parser.parse_args()
args.column = args.column - 1

# global variables
dirs = os.listdir(args.directory)
dictionary = {}
columnNames = []

for fileName in dirs:
    if fileName != args.out:
        if args.name in fileName:

            # Get Sample Name
            fileNameSplit = fileName.split(".")
            sampleName = fileNameSplit[0]

            print("Processing "+fileName, file = sys.stderr)

            columnNames.append(sampleName)

            # Process file
            with open(args.directory+"/"+fileName) as f:
                for line in f:
                    split = line.strip().split("\t")

                    rowName = split[0]

                    value = int(split[args.column])
                    if value > args.minimum:
                        if rowName in dictionary:
                            dictionary[rowName][sampleName] = value
                        else:
                            dictionary[rowName] = {sampleName : value}

print("Printing results", file = sys.stderr)
outFile = open(args.out, 'w')

for sampleName in columnNames:
    outFile.write(sampleName+"\t")
outFile.write("\n")

for rowName in sorted(dictionary.keys()):
    outFile.write(rowName)
    for sampleName in columnNames:
        if sampleName in dictionary[rowName]:
            outFile.write("\t"+str(dictionary[rowName][sampleName]))
        else:
            outFile.write("\t0")
    outFile.write("\n")
