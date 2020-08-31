import os
import sys
import argparse
from Bio import SeqIO

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Reciprocal best blast, yo.')

parser.add_argument('-a', '--proteinA',
                    help="First protein fasta file",
                    required=True)

parser.add_argument('-b', '--proteinB',
                    help="Second protein fasta file",
                    required=True)

parser.add_argument('-l', '--length',
                    help="Minimum homology length",
                    default=0.5)

parser.add_argument('-e', '--eval',
                    help="Minimum e-value",
                    default=1e-15)

parser.add_argument('-i', '--percentID',
                    help="Minimum percent identity",
                    default=0.5)

args = parser.parse_args()
args.length = float(args.length)
args.eval = float(args.eval)
args.percentID = float(args.percentID)
if args.percentID < 1:  # blast given as percentage rather than ratio
    args.percentID = args.percentID * 100

########### Process File A ###########

os.system('makeblastdb -in '+args.proteinA+' -dbtype prot')
aDict = {}
for seq_record in SeqIO.parse(args.proteinA, "fasta"):
    seq_record.seq = seq_record.seq.rstrip("*")
    aDict[seq_record.id] = {'length': len(
        seq_record.seq), "hit": "", "percentLength": 0, "eval": 1, "percentID": 0, "reciprocal": 0}

########### Process File B ###########

os.system('makeblastdb -in '+args.proteinB+' -dbtype prot')
bDict = {}
for seq_record in SeqIO.parse(args.proteinB, "fasta"):
    seq_record.seq = seq_record.seq.rstrip("*")
    bDict[seq_record.id] = {'length': len(
        seq_record.seq), "hit": "", "percentLength": 0, "eval": 1, "percentID": 0, "reciprocal": 0}

########### RUN BLASTs ###########

os.system('blastp -query '+args.proteinA+' -db ' +
          args.proteinB+' -evalue 1e-7 -outfmt 6 -out a_against_b.blastp')
os.system('blastp -query '+args.proteinB+' -db ' +
          args.proteinA+' -evalue 1e-7 -outfmt 6 -out b_against_a.blastp')

########### PROCESS BLASTs ###########

for line in open('a_against_b.blastp', 'rt'):
    line = line.rstrip()
    [query, subject, percent_id, alignment_length, mismatches, gap_openings, query_start,
        query_end, subject_start, subject_end, E_value, bit_score] = line.split("\t")

    # add hits satisfying criteria to dictionary
    percentLength = int(alignment_length) / aDict[query]['length']
    if percentLength >= args.length:
        if float(E_value) <= args.eval:
            if float(percent_id) >= args.percentID:
                if percentLength >= aDict[query]['percentLength']:
                    if float(E_value) <= aDict[query]['eval']:
                        aDict[query]['hit'] = subject
                        aDict[query]['percentLength'] = percentLength
                        aDict[query]['alignLength'] = int(alignment_length)
                        aDict[query]['eval'] = float(E_value)
                        aDict[query]['percentID'] = float(percent_id)

for line in open('b_against_a.blastp', 'rt'):
    line = line.rstrip()
    [query, subject, percent_id, alignment_length, mismatches, gap_openings, query_start,
        query_end, subject_start, subject_end, E_value, bit_score] = line.split("\t")

    # add hits satisfying criteria to dictionary
    percentLength = int(alignment_length) / bDict[query]['length']
    if percentLength >= args.length:
        if float(E_value) <= args.eval:
            if float(percent_id) >= args.percentID:
                if percentLength >= bDict[query]['percentLength']:
                    if float(E_value) <= bDict[query]['eval']:
                        bDict[query]['hit'] = subject
                        bDict[query]['percentLength'] = percentLength
                        bDict[query]['alignLength'] = int(alignment_length)
                        bDict[query]['eval'] = float(E_value)
                        bDict[query]['percentID'] = float(percent_id)

########### FIND RECIPROCAL HITS ###########

for query in aDict:
    subject = aDict[query]['hit']
    if subject != "":
        if query == bDict[subject]['hit']:
            aDict[query]['reciprocal'] = 1
            bDict[subject]['reciprocal'] = 1

for query in sorted(aDict.keys()):
    if aDict[query]['reciprocal'] == 1:
        subject = aDict[query]['hit']
        print(query, aDict[query]['length'], aDict[query]['alignLength'],
              aDict[query]['percentID'], aDict[query]['eval'], sep="\t", end="\t")
        print(subject, bDict[subject]['length'], bDict[subject]['alignLength'],
              bDict[subject]['percentID'], bDict[subject]['eval'], sep="\t")
    else:
        print(query, "-", "-", "-", "-", "None", "-", "-", "-", "-", sep="\t")

for query in sorted(bDict.keys()):
    if bDict[query]['reciprocal'] != 1:
        print("None", "-", "-", "-", "-", query, "-", "-", "-", "-", sep="\t")
