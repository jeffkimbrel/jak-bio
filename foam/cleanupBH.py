import sys

lines = [line.strip() for line in open(sys.argv[1])]

for line in list(lines):
    split = line.split()
    if len(split) > 0:
        underscoreSplit = split[3].split("_")
        split[3] = underscoreSplit[0]
        print("\t".join(split))