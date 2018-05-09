import sys
from Bio import SeqIO

frames = [1,2,3]

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    for frame in frames:
        absFrame = frame
        if frame > 0:
            seq = seq_record.seq
        elif frame < 0:
            seq = seq_record.seq.reverse_complement()
            absFrame = abs(frame)
        
        # trim seq to modulus 3 to remove the biopython error
        trim = len(seq[absFrame-1:]) % 3
        trim = 0 - trim
        if trim != 0:
            seq = seq[:trim]
        
        print(">"+seq_record.id+"_"+str(frame))      
        print(seq[absFrame-1:].translate())
