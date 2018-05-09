import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'Takes a list of genes and cazy annotations, and produces summary statistics for classes and berlemont families')

parser.add_argument('-f', '--file', help = "List of genes", required = True)
parser.add_argument('-t', '--type', default = "f", help = "(f)amily, (c)lass or (b)erlemont" )
parser.add_argument('-c', '--column', default = "2", help = "Column number with CAZy annotations" )
parser.add_argument('--markdown', '-m', action='store_true', help='Print as markdown' )

args = parser.parse_args()
args.type = args.type.lower()
args.column = int(args.column) - 1

fileTypes = {"f" : "Family", "c" : "Class" , "b" : "Berlemont"}

lines = [line.strip() for line in open(args.file)]

dictionary = {}

berlemontFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/bio/cazy/berlemont_types.txt")]
berlemont = {}
for line in berlemontFile:
    line = line.strip()
    split = line.split("\t")

    berlemont[split[0]] = split[1]

for line in lines:
    cazy = line.split("\t")[args.column]

    if args.type == "f":
        dictionary[cazy] = dictionary.get(cazy, 0) + 1

    if args.type == "c":
        cazy = cazy.rstrip('1234567890')
        dictionary[cazy] = dictionary.get(cazy, 0) + 1

    if args.type == "b":
        if cazy in berlemont:
            cazy = berlemont[cazy]
            dictionary[cazy] = dictionary.get(cazy, 0) + 1
        else:
            dictionary["NA"] = dictionary.get("NA", 0) + 1

separator = "\t"
if args.markdown == True:
    separator = " | "
    print(fileTypes[args.type], "|", "Count")
    print("--- | --- ")

for cazy in sorted(dictionary.keys()):
    print(cazy, dictionary[cazy], sep = separator)
