import sys
from Bio import SeqIO

# this version appends abundance information to the end of the fasta header

class cluster:
    clusterList = []

    def __init__(self, clusterNumber, rep, members,abundance):
        cluster.clusterList.append(self)
        self.cluster = clusterNumber
        self.rep = rep
        self.members = members
        self.abundance = abundance

def processHeader(header):
    repSplit = header.split(" ")
    currentHeader = repSplit[1][1:-3]

    abundance = 1

    split = currentHeader.split(";")
    if len(split) > 1:
        split2 = split[1].split("=")
        abundance = int(split2[1])
        #currentHeader = split[0]

    return(currentHeader,abundance)

# first, read in the cluster file and get the abundances of each

file = [line.strip() for line in open(sys.argv[1])]

currentCluster = "none"
currentRep = ""
currentMembers = []
currentAbundance = 0
largestCluster = 0

for line in file:
    if line.startswith(">"):
        if currentCluster != "none":
            cluster(currentCluster,currentRep,currentMembers,currentAbundance)

        clusterSplit = line.split(" ")
        currentAbundance = 0
        currentCluster = clusterSplit[1]
    elif line.startswith("0"):
        currentRep,abundance = processHeader(line)
        currentMembers = [currentRep]
        currentAbundance = currentAbundance + abundance
    else:
        currentMember,abundance = processHeader(line)
        currentMembers.append(currentMember)
        currentAbundance = currentAbundance + abundance

cluster(currentCluster,currentRep,currentMembers,currentAbundance)

# get largest cluster
for clusterNumber in cluster.clusterList:
    if clusterNumber.abundance > largestCluster:
        largestCluster = clusterNumber.abundance

# fasta

currentLength  = largestCluster
totalAbundance = 0

#record = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
record = SeqIO.index(sys.argv[2], "fasta")

while currentLength > 0:
    for clusterNumber in cluster.clusterList:
        if clusterNumber.abundance == currentLength:

            totalAbundance = totalAbundance + currentLength

            seq = record[clusterNumber.rep].seq

            # remove old abundance tag
            newClusterName = clusterNumber.rep.split(";")

            #print(">cd" + str(clusterNumber.cluster) + " a=" + str(clusterNumber.abundance) + "\n" + seq)
            print(">cd" + str(clusterNumber.cluster) + "\n" + seq)
            #print(record[clusterNumber.rep].format("fasta"))

    currentLength = currentLength - 1
