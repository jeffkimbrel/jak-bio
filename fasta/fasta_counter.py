from Bio import SeqIO
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = 'Given a fasta file, returns the total or subsequence nt/aa count.')
parser.add_argument('-f', '--fasta', help = "fasta file", required = True)
parser.add_argument('--total', '-t', action = 'store_true', help = 'Print total of all subsequences' )
parser.add_argument('--sequence', '-s', action = 'store_true', help = 'Print total for each subsequence' )
args = parser.parse_args()

total = 0

## FUNCTIONS ###################################################################

def addToTotal(count):
    global total
    total += count

def printSequenceTotal(id, count):
    if args.sequence == True:
        print(id, count, sep = "\t")

def printTotal():
    if args.total == True:
        print("Total", total, sep = "\t")

def getInfo(seq_record):
    id = str(seq_record.id)
    count = len(str(seq_record.seq))
    return(id, count)

## ANALYZE #####################################################################

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    id, count = getInfo(seq_record)
    addToTotal(count)
    printSequenceTotal(id, count)

printTotal()

if args.total == False:
    if args.sequence == False:
        print("\nYou should probably use the -s or -t flags to see some results...\n")
