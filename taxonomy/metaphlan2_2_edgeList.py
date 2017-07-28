import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-f', '--file',
    help = "MetaPhlan2 output",
    required = True)

parser.add_argument('-k', '--kingdom',
    help = "Kingdoms to include (ABEV)",
    default = "AB")

parser.add_argument('-c', '--cutoff',
    default = 0.01,
    help = "Percent cutoff of taxa to include" )

args = parser.parse_args()
args.kingdom = args.kingdom.upper()

args.cutoff = float(args.cutoff)
if float(args.cutoff) >= 1:
    args.cutoff /= 100

# Functions

def includeTaxaLine(line):
    value = False

    if "A" in args.kingdom:
        if "k__Archaea" in line:
            value = True
    if "B" in args.kingdom:
        if "k__Bacteria" in line:
            value = True
    if "E" in args.kingdom:
        if "k__Eukaryota" in line:
            value = True
    if "V" in args.kingdom:
        if "k__Viruses" in line:
            value = True

    return(value)


# Process data

f = open(args.file,'r')
lines = f.readlines()
f.close()

for line in lines:
    if not line.lstrip().startswith('#'):
        line = line.strip()
        taxa,score = line.split("\t")
        score = float(score) / 100

        if includeTaxaLine(line):
            if score >= args.cutoff:
                taxaList = taxa.split("|")
                taxaList.insert(0, "Root")

                print(*taxaList[-2:], sep = "\t", end = "\t")
                print(score)
