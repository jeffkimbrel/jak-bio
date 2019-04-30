import sys
import os
import argparse
from Bio import SeqIO

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-a', '--abundance',
    required = True,
    help = "File with raw counts" )

parser.add_argument('-f', '--fasta',
    required = True,
    help = "Fasta files for length determination" )

parser.add_argument('-c', '--column',
    default = 2,
    type = int,
    help = "Column with raw counts (1-based)" )

parser.add_argument('-g', '--geneNames',
    default = 1,
    help = "Column with gene names (1-based)" )

args = parser.parse_args()
args.column = args.column - 1
args.geneNames = args.geneNames - 1

## FUNCTIONS ###################################################################

def cleanGeneName(gene):
    split = gene.split(" ")
    return(split[0])

## FASTA #######################################################################

genes = {}

for seq_record in SeqIO.parse(args.fasta, "fasta"):
    genes[seq_record.id] = {'LENGTH' : len(str(seq_record.seq)) / 1000, 'RAW' : 0, "TPM" : 0, 'RPK' : 0}

## ABUNDANCE ###################################################################

totalRPK = 0

abundanceFile = [line.strip() for line in open(args.abundance)]

for line in abundanceFile:
    split = line.split("\t")

    gene = cleanGeneName(split[args.geneNames])
    raw = split[args.column]

    if raw.isdigit():
        raw = float(raw)
        genes[gene]['RAW'] = raw
        genes[gene]['RPK'] = raw / genes[gene]['LENGTH']
        totalRPK += genes[gene]['RPK']

perMillion = totalRPK / 1000000

## CALCULATE TPM ###############################################################

for gene in genes:
    genes[gene]['TPM'] = genes[gene]['RPK'] / perMillion

    #print(gene, genes[gene]['LENGTH'], genes[gene]['RAW'], genes[gene]['RPK'], genes[gene]['TPM'], sep = "\t")
    print(gene, genes[gene]['TPM'], sep = "\t")
