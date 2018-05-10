import sys
import argparse

parser = argparse.ArgumentParser(description = 'Takes a cd-hit cluster file (.clstr) and a FOAM.BH file, and assigns the FOAM annotations to the cluster members.')

parser.add_argument('-f', '--foam',
    help = "FOAM.BH file",
    required = True)

parser.add_argument('-c', '--clstr',
    help = "cd-hit .clstr file",
    required = True)

parser.add_argument('--split', '-s',
    action = 'store_true',
    help = 'Split multi-KEGGs to their own line' )

parser.add_argument('--ontology', '-o',
    action = 'store_true',
    help = 'Include complete ontology for each KO. Overrides -s option.' )

args = parser.parse_args()

##### CLASSES/FUNCTIONS ########################################################

class cluster:
    clusterList = []

    def __init__(self, clusterNumber, rep, members):
        cluster.clusterList.append(self)
        self.cluster = clusterNumber
        self.rep = rep
        self.members = members

def processHeader(header):
    repSplit = header.split(" ")
    currentHeader = repSplit[1][1:-3]
    split = currentHeader.split(";")
    if len(split) > 1:
        split2 = split[1].split("=")

    return(currentHeader)

##### FOAM FILE ################################################################
foam = {}
foamFile = [line.strip() for line in open(args.foam)]

for line in foamFile:
    split = line.split("\t")
    foam[split[0]] = split[3].replace("KO:", "")

##### ONTOLOGY FILE ############################################################
ontology = {}
ontologyFile = [line.strip() for line in open("/Users/kimbrel1/Dropbox/scripts/FOAM/FOAM-onto_rel1.tsv")]

for line in ontologyFile:
    split = line.split("\t")

    if split[4] in ontology:
        ontology[split[4]].append(split[0:4])
    else:
        ontology[split[4]] = [split[0:4]]

##### CLSTR FILE ###############################################################

clstr = [line.strip() for line in open(args.clstr)]

currentCluster = "none"
currentRep = ""
currentMembers = []

for line in clstr:
    if line.startswith(">"):
        if currentCluster != "none":
            cluster(currentCluster, currentRep, currentMembers)

        clusterSplit = line.split(" ")
        currentAbundance = 0
        currentCluster = clusterSplit[1]
    elif '*' in line:
        currentRep = processHeader(line)
        currentMembers = [currentRep]
    else:
        currentMember = processHeader(line)
        currentMembers.append(currentMember)

# catches the last one
cluster(currentCluster, currentRep, currentMembers)


print("ORF", "FOAM", "SCALE", "CLUSTER", "L1", "L2", "L3", "L4", sep = "\t")

for clusterName in cluster.clusterList:
    if clusterName.rep in foam:

        for member in clusterName.members:
            if args.ontology == True:
                split = foam[clusterName.rep].split(",")
                for ko in split:

                    scale = 1 / len(ontology[ko])

                    for koOntology in ontology[ko]:
                        print(member, ko, scale, clusterName.cluster, *koOntology, sep = "\t")
            elif args.split == True:
                split = foam[clusterName.rep].split(",")

                scale = 1 / len(split)

                for ko in split:
                    print(member, ko, scale, clusterName.cluster, "", "", "", "", sep = "\t")
            else:
                print(member, foam[clusterName.rep], 1.0, clusterName.cluster, "", "", "", "", sep = "\t")
    else:
        for member in clusterName.members:
            print(member, "NONE", 1.0, clusterName.cluster, "", "", "", "", sep = "\t")
