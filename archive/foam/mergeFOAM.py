import glob
import sys

S = {} # key = sample, value = KEGG:counts
K = {} # key = kegg, value = kegg
level = "KO" # L1,L2,L3,L4,KO

# ONTOLOGY
ontology = {}
notInOntology = {}

ontFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/FOAM/FOAM-onto_rel1.tsv")]

for line in ontFile:
    L1,L2,L3,L4,KO = line.split("\t")

    if KO in ontFile:
        ontology[KO]['L1'].append[L1]
        list(set(ontology[KO]['L1']))
        ontology[KO]['L2'].append[L2]
        list(set(ontology[KO]['L2']))
        ontology[KO]['L3'].append[L3]
        list(set(ontology[KO]['L3']))
        ontology[KO]['L4'].append[L4]
        list(set(ontology[KO]['L4']))
        ontology[KO]['KO'].append[KO]
        list(set(ontology[KO]['KO']))

    else:
        ontology[KO] = {'L1' : L1, 'L2' : L2, 'L3' : L3, 'L4' : L4, 'KO' : KO}

for f in glob.glob('*.ko.txt'):

    S[f] = {}
    lines = [line.strip() for line in open(f)]

    for line in lines:
        kegg, count = line.split()

        if kegg not in ontology:
            notInOntology[kegg] = kegg

        if level == 'KO':
            K[kegg] = kegg
            S[f][kegg] = count
        #elif level == 'L1':


sortedK = sorted(list(K))
sortedS = sorted(list(S))


## HEADER
print(level, end="")
for sample in sortedS:
    split = sample.split(".")
    print("\t"+split[0],end="")
    #print("\t"+sample, end = "")

print()

## COUNTS
for kegg in sortedK:
    print(kegg, end="")
    for sample in sortedS:
        #print(S[sample][kegg])
        print("\t"+str(S[sample].get(kegg, 0)),end="")
    print()


# NON-ONTOLOGY
for kegg in notInOntology:
    print(kegg)
