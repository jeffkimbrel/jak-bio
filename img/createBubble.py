import glob
import sys

phylo = {}
depth = {}
lengthGC = {}

level = 2

for f in glob.glob('*.txt'):
    if "phylodist" in f:
        lines = [line.strip() for line in open(f)]

        for line in lines:
            split = line.split()

            taxa = split[4].split(";")
            phylo[split[0]] = taxa[level]
            print(taxa[level])
