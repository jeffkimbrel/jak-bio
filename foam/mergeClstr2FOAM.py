import sys
import argparse

parser = argparse.ArgumentParser(description = 'X')

parser.add_argument('-f', '--foam',
    help = "FOAM.BH file",
    required = True)

parser.add_argument('-c', '--clstr',
    help = "cd-hit .clstr file",
    required = True)

parser.add_argument('--split', '-s',
    action = 'store_true',
    help = 'Split multi-KEGGs to their own line' )

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

for clusterName in cluster.clusterList:
    if clusterName.rep in foam:

        for member in clusterName.members:
            if args.split == True:
                split = foam[clusterName.rep].split(",")
                for ko in split:
                    print(member, ko, clusterName.cluster, sep = "\t")
            else:
                print(member, foam[clusterName.rep], clusterName.cluster, sep = "\t")
    else:
        for member in clusterName.members:
            print(member, "NONE", clusterName.cluster, sep = "\t")
