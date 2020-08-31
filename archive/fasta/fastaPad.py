import sys
from Bio import SeqIO

longest = 0

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    if len(seq_record.seq) > longest:
        longest = len(seq_record.seq)

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    seq_record.seq = seq_record.seq + longest*"X"
    print(">"+seq_record.id+"\n"+seq_record.seq[0:longest])