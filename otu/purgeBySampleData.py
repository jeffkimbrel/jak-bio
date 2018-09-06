import os
import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'Howdy!', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-s', '--sampledata',
    required = True,
    help = "Sample Data" )

parser.add_argument('-t', '--table',
    required = True,
    help = "ASV (OTU) table" )

parser.add_argument('-ms', '--minimumSamples',
    default = 2,
    type = int,
    help = "Must be found in at least this many samples within a group" )

parser.add_argument('-mr', '--minimumReads',
    default = 1,
    type = int,
    help = "Minimum number of reads to be considered 'found'" )

parser.add_argument('-c', '--columns',
    required = True,
    help = "Column numbers (can be multiple, space-separated) to group samples by",
    nargs = "*",
    type = int)

args = parser.parse_args()

## Functions ###################################################################

def getGroupName(line):
    split = line.split("\t")
    groupName = []
    for c in args.columns:
        groupName.append(split[c])
    return("_".join(groupName))

## Classes #####################################################################

groups = {}
class Group:
    def __init__(self, groupName, sample):
        self.groupName = groupName
        self.samples = [sample]
        self.count = 1

    def addSample(self, sample):
        self.samples.append(sample)
        self.count += 1

samples = {}
class Sample:
    def __init__(self, groupName, sample):
        self.sample = sample
        self.groupName = groupName
        self.asv = {}

    def addCount(self, ASV, count):
        self.asv[ASV] = count

asvs = {}
class ASV:
    def __init__(self, asv):
        self.asv = asv
        self.counts = {}

    def addCount(self, passedSamples, sample, count):
        if passedSamples >= args.minimumSamples:
            self.counts[sample] = count
        else:
            self.counts[sample] = 0

    def printASV(self):
        asvSum = 0
        counter = 1
        while counter < len(ASVheader):
            sample = ASVheader[counter]
            asvSum += self.counts[sample]

            counter += 1
        if asvSum > 0:
            print(self.asv, end = "\t")
            counter = 1
            while counter < len(ASVheader):
                sample = ASVheader[counter]
                print(self.counts[sample], end = "\t")
                counter += 1
            print()




## Group Samples ###############################################################

sampleDataFile = [line.strip() for line in open(args.sampledata)]
sampleDataFile.pop(0) # remove header

for line in sampleDataFile:
    split = line.split("\t")

    groupName = getGroupName(line)

    ## create Sample
    samples[split[0]] = Sample(groupName, split[0])

    ## create or add to a Group
    if groupName in groups:
        groups[groupName].addSample(split[0])
    else:
        groups[groupName] = Group(groupName, split[0])

# print summary
print("GROUP", "N", "SAMPLE", sep = "\t", file = sys.stderr)
for group in groups:
    print(groups[group].groupName, groups[group].count, groups[group].samples, sep = "\t", file = sys.stderr)
print(file = sys.stderr)

## Process ASVs ################################################################

tableFile = [line.strip() for line in open(args.table)]
ASVheader = tableFile.pop(0).split("\t")

for line in tableFile:
    split = line.split("\t")
    asv = split[0]

    asvs[asv] = ASV(asv)

    counter = 1
    while counter < len(split):
        sample = ASVheader[counter]
        count = int(split[counter])

        if count < args.minimumReads: # Zero out low counts
            count = 0

        samples[sample].addCount(asv, count)
        counter += 1

## Find Supported Groups #######################################################

for group in groups:

    for asv in asvs:
        passedSamples = 0
        for sample in groups[group].samples:
            sampleCount = samples[sample].asv[asv]

            if sampleCount > 0: #already set minimum reads to zero above
                passedSamples += 1

        # now go again based on the passedSample Count
        for sample in groups[group].samples:
            sampleCount = samples[sample].asv[asv]

            asvs[asv].addCount(passedSamples, sample, sampleCount)

## Print New Table #############################################################

counter = 1
while counter < len(ASVheader):
    print(ASVheader[counter], end = "\t")
    counter += 1
print()

for asv in sorted(asvs.keys()):
    asvs[asv].printASV()
