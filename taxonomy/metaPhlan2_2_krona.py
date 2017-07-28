import sys

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

for line in lines:
    if not line.lstrip().startswith('#'):
        line = line.strip()
        taxa, score = line.split("\t")

        taxaList = taxa.split("|")

        if len(taxaList) == 7:

            print(score, *taxaList, sep = "\t")
