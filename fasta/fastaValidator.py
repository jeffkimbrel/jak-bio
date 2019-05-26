import sys
from Bio import SeqIO

fileName = sys.argv[1]
handle = open(fileName, "r")

seqCount = 0
emptyHeaders = 0
dupHeaders = {}

for record in SeqIO.parse(handle, "fasta"):
    seqCount += 1
    if record.id == "":
        emptyHeaders += 1
    else:
        dupHeaders[record.id] = dupHeaders.get(record.id, 0) + 1

dupHeaderFound = 0
print("Duplicate Headers: ")

for header in dupHeaders:
    if dupHeaders[header] > 1:
        print(header, dupHeaders[header], sep="\t")
        dupHeaderFound = 1

if dupHeaderFound == 0:
    print("None")

print("\nSeqs: ",seqCount,"\nEmpty Headers: ", emptyHeaders)
