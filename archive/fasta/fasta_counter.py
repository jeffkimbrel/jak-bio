from Bio import SeqIO
import argparse

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Given a fasta file, \
    returns the total or subsequence nt/aa count.')

parser.add_argument('-f',
                    '--fasta',
                    nargs='*',
                    help="fasta file",
                    required=True)
parser.add_argument('--total',
                    '-t',
                    action='store_true',
                    help='Print total of all subsequences')
parser.add_argument('--sequence',
                    '-s',
                    action='store_true',
                    help='Print total for each subsequence')

args = parser.parse_args()

# FUNCTIONS ###################################################################


def printSequenceTotal(id, count):
    if args.sequence is True:
        print(id, count, sep="\t")


def getInfo(seq_record):
    id = str(seq_record.id)
    count = len(str(seq_record.seq))
    return(id, count)

# ANALYZE #####################################################################


for fasta_file in args.fasta:
    total = 0
    for seq_record in SeqIO.parse(fasta_file, "fasta"):

        id, count = getInfo(seq_record)
        total += count
        printSequenceTotal(id, count)

    if args.total is True:
        print(fasta_file, total, sep="\t")

if args.total is False:
    if args.sequence is False:
        print("\nYou should probably use the -s or -t flags to see some results...\n")
