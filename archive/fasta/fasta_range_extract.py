import sys
from Bio import SeqIO
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Extract a sequence range from a fasta file')

parser.add_argument('-f', '--fasta',
    help="fasta file",
    required=True)

parser.add_argument('-b', '--begin',
    help="Beginning nucleotide position (1-based)",
    type = int,
    required=True)

parser.add_argument('-s', '--sequence',
    default = "",
    help="Only extract from sequence with this header")

parser.add_argument('-e', '--end',
    help="Ending nucleotide position (1-based)",
    type = int,
    required=True)

parser.add_argument('--revcom', '-rc',
    action = 'store_true',
    help = 'Reverse Complement?' )

args = parser.parse_args()

args.begin -= 1

## MAIN ########################################################################

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    id = str(seq_record.id)

    if len(args.sequence) > 0:
        if id == args.sequence:
            id = id + "_" + str(args.begin + 1) + "-" + str(args.end)

            if args.revcom == True:
                print(">" + id + "rc")
                print(seq_record.seq[args.begin:args.end].reverse_complement())
            else:
                print(">" + id)
                print(seq_record.seq[args.begin:args.end])
    else:
        id = id + "_" + str(args.begin + 1) + "-" + str(args.end)

        if args.revcom == True:
            print(">" + id + "rc")
            print(seq_record.seq[args.begin:args.end].reverse_complement())
        else:
            print(">" + id)
            print(seq_record.seq[args.begin:args.end])
