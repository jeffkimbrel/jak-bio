import sys
from Bio import SeqIO
import argparse
from multiprocessing import Pool

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-f', '--fasta',
                    required=True)

parser.add_argument('--frames',
                    nargs="*",
                    default=[1])

args = parser.parse_args()

###

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    for frame in args.frames:
        frame = int(frame)
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
        print(seq[absFrame-1:].translate(table=11))
