import sys
from Bio import SeqIO
import re

fileName = sys.argv[1]

handle = open(fileName, "rU")

for record in SeqIO.parse(handle, "fasta"):
    seqOrig = str(record.seq)
    seqNew = re.sub("[^a-zA-Z]+","",seqOrig)
    
    print(">"+str(record.id)+"\n"+str(seqNew))
    