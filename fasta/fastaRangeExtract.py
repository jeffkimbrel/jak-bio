import sys
from Bio import SeqIO

start = int(sys.argv[2]) - 1
stop = int(sys.argv[3])

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    print(">"+seq_record.id)
    print(seq_record.seq[start:stop])
