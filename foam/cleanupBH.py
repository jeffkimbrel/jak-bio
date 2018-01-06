# Jeff Kimbrel, jakpot@gmail.com

import sys


orphanFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/FOAM/orphanKOs.txt")]

orphans = {}

for line in orphanFile:
    orphans[line] = {}




lines = [line.strip() for line in open(sys.argv[1])]

for line in list(lines):
    split = line.split()
    if len(split) > 0:
        underscoreSplit = split[3].split("_")

        # remove orphans
        nonOrphanList = []
        commaSplit = underscoreSplit[0].split(",")
        for ko in commaSplit:

            kegg = ""

            if "KO:" in ko:
                colonSplit = ko.split(":")
                kegg = colonSplit[1]
            else:
                kegg = ko

            if kegg not in orphans:
                nonOrphanList.append("KO:"+kegg)

        nonOrphans = ",".join(nonOrphanList)

        #split[3] = underscoreSplit[0]
        split[3] = nonOrphans

        print("\t".join(split))
