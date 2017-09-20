import sys

levelSelect = 1 # 1=L1, 2=L2, etc.
levelSelect = levelSelect - 1
ontology = {}
ontFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/FOAM/FOAM-onto_rel1.tsv")]
junk = ontFile.pop(0)

newTable = {}

for line in ontFile:
    L1,L2,L3,L4,KO = line.split("\t")

    L1 = "L1-"+L1

    if L2 == "":
        L2 = L1
    else:
        L2 = L1+"_L2-"+L2

    if L3 == "":
        L3 = L2
    else:
        L3 = L2+"_L3-"+L3

    if L4 == "":
        L4 = L3
    else:
        L4 = L3+"_L4-"+L4

    lineOnt = (L1,L2,L3,L4)

    if KO in ontology:
        ontology[KO].append(lineOnt)
    else:
        ontology[KO] = []
        ontology[KO].append(lineOnt)

lines = [line.strip() for line in open(sys.argv[1])]
#print(lines.pop(0))

for line in list(lines):

    split = line.split()
    KO = split.pop(0)

    levelList = []
    if KO in ontology: # ignores KOs not found in the ontology
        for level in ontology[KO]:
            levelList.append(level[levelSelect])

    levelList = list(set(levelList))
    for level in levelList:
        if level in newTable:
            count = 0

            while count < len(split):
                newTable[level][count] = "{0:.2f}".format(float(newTable[level][count]) + float(split[count]))
                count = count + 1
        else:
            newTable[level] = []
            count = 0
            while count < len(split):
                newTable[level].append(split[count])
                count = count + 1

for level in newTable:
    print(level,*newTable[level],sep="\t")
