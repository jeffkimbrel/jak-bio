import sys
from Bio import SeqIO

class cluster:
    clusterList = []
    
    def __init__(self, clusterNumber, rep, members):
        cluster.clusterList.append(self)
        self.cluster = clusterNumber
        self.rep = rep
        self.members = members
        self.abundance = len(members)  
        
# first, read in the cluster file and get the abundances of each

file = [line.strip() for line in open(sys.argv[1])]

currentCluster = "none"
currentRep = ""
currentMembers = []
largestCluster = 0

for line in file:
    if line.startswith(">"):
        
        if currentCluster != "none":
            cluster(currentCluster,currentRep,currentMembers)
                  
        clusterSplit = line.split(" ")
        currentCluster = clusterSplit[1]
    elif line.startswith("0"):
        repSplit = line.split(" ")
        currentRep = repSplit[1][1:-3]
        currentMembers = [currentRep]
    else:
        memberSplit = line.split(" ")
        currentMember = memberSplit[1][1:-3]
        currentMembers.append(currentMember)
        
cluster(currentCluster,currentRep,currentMembers)
        
# get largest cluster        
for clusterNumber in cluster.clusterList:
    if clusterNumber.abundance > largestCluster:
        largestCluster = clusterNumber.abundance

# fasta

currentLength  = largestCluster

record = SeqIO.index(sys.argv[2], "fasta")

while currentLength > 0:
    for clusterNumber in cluster.clusterList:
        if clusterNumber.abundance == currentLength:
            print(record[clusterNumber.rep].format("fasta"))
                             
    currentLength = currentLength - 1
    