import os
import sys
import argparse
from Bio import SeqIO
from natsort import natsorted
import pandas as pd

from jakomics import utilities, blast, colors

import jak_utils

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Reciprocal best blast, yo.')

parser.add_argument('-a', '--proteinA',
                    help="First protein fasta file",
                    required=True)

parser.add_argument('-b', '--proteinB',
                    help="Second protein fasta file",
                    required=True)

parser.add_argument('-o', '--out',
                    help="Out file",
                    type=str,
                    default="reciprocal_blast.txt")

parser.add_argument('-e', '--eval',
                    help="Minimum e-value, set at the actual blast step",
                    type=float,
                    default=1e-7)

parser.add_argument('-l', '--length',
                    help="Minimum homology length",
                    type=float,
                    default=0.0)

parser.add_argument('-i', '--percentID',
                    help="Minimum percent identity",
                    type=float,
                    default=0.0)

parser.add_argument('-t', '--threads',
                    help="Blast threads",
                    type=int,
                    default=1)

args = parser.parse_args()

if args.percentID < 1.0:  # blast given as ratio rather than percentage
    args.percentID = args.percentID * 100.0

# FUNCTIONS ####################################################################


def make_gene_dict(A, B):
    gene_dict = {}
    gene_dict['A'] = SeqIO.to_dict(SeqIO.parse(A, "fasta"))
    gene_dict['B'] = SeqIO.to_dict(SeqIO.parse(B, "fasta"))

    return gene_dict


def get_gene_length(gene_set, locus_tag):
    return len(gene_dict[gene_set][locus_tag].seq)


# MAIN #########################################################################


if __name__ == "__main__":
    jak_utils.header()
    print(f"{colors.bcolors.RED}Code is under development{colors.bcolors.END}")

    # collect gene information
    gene_dict = make_gene_dict(args.proteinA, args.proteinB)
    best_dict = {'A': {},
                 'B': {}
                 }

    # Do BLASTS
    a_against_b = blast.run_blast(type="prot",
                                  q=args.proteinA,
                                  db=args.proteinB,
                                  e=args.eval,
                                  threads=args.threads,
                                  echo=True,
                                  make=True)

    b_against_a = blast.run_blast(type="prot",
                                  q=args.proteinB,
                                  db=args.proteinA,
                                  e=args.eval,
                                  threads=args.threads,
                                  echo=True,
                                  make=True)

    # process A results
    for locus_tag, gene_info in gene_dict['A'].items():
        if locus_tag in a_against_b:
            best = a_against_b[locus_tag][0]  # first hit will be the "best" hit
            if best.alignment_length / len(gene_info.seq) >= args.length:
                if best.percent >= args.percentID:
                    best_dict['A'][locus_tag] = {'locus_tag': locus_tag,
                                                 'subject': best.subject,
                                                 'eval': best.eval,
                                                 'alignment_length': best.alignment_length,
                                                 'mismatches': best.mismatches,
                                                 'bit_score': best.bit_score,
                                                 'query_start': best.query_start,
                                                 'query_end': best.query_end,
                                                 'subject_start': best.subject_start,
                                                 'subject_end': best.subject_end,
                                                 'already_paired': False
                                                 }

    # process B results
    for locus_tag, gene_info in gene_dict['B'].items():
        if locus_tag in b_against_a:
            best = b_against_a[locus_tag][0]  # first hit will be the "best" hit
            if best.alignment_length / len(gene_info.seq) >= args.length:
                if best.percent >= args.percentID:
                    best_dict['B'][locus_tag] = {'locus_tag': locus_tag,
                                                 'subject': best.subject,
                                                 'eval': best.eval,
                                                 'alignment_length': best.alignment_length,
                                                 'mismatches': best.mismatches,
                                                 'bit_score': best.bit_score,
                                                 'query_start': best.query_start,
                                                 'query_end': best.query_end,
                                                 'subject_start': best.subject_start,
                                                 'subject_end': best.subject_end,
                                                 'already_paired': False
                                                 }

    df = pd.DataFrame(columns=['RECIPROCAL', 'A_LOCUS', 'B_LOCUS', 'A_LENGTH', 'A_COORDS', 'A_EVAL', 'A_ALIGNMENT_LENGTH',
                               'A_MISMATCHES', 'A_BITSCORE', 'B_LENGTH', 'B_COORDS', 'B_EVAL', 'B_ALIGNMENT_LENGTH', 'B_MISMATCHES', 'B_BITSCORE'])

    # process A hits
    for locus_tag, a_hit in best_dict['A'].items():
        # print(f"---{locus_tag}-==")
        if a_hit['already_paired'] == False:
            if a_hit['subject'] in best_dict['B']:
                b_hit = best_dict['B'][a_hit['subject']]

                if b_hit['subject'] == locus_tag:
                    if b_hit['already_paired'] == False:
                        best_dict['A'][locus_tag]['already_paired'] = True
                        best_dict['B'][a_hit['subject']]['already_paired'] = True
                        # print(a_hit)
                        # print(b_hit)
                        #
                        # print(f'{colors.bcolors.GREEN}Besties!!{colors.bcolors.END}')

                        df = df.append(
                            pd.Series(data={
                                'RECIPROCAL': "AB",
                                'A_LOCUS': locus_tag,
                                'B_LOCUS': a_hit['subject'],
                                'A_LENGTH': get_gene_length('A', locus_tag),
                                'A_COORDS': f"{a_hit['query_start']}-{a_hit['query_end']}",
                                'A_EVAL': a_hit['eval'],
                                'A_ALIGNMENT_LENGTH': a_hit['alignment_length'],
                                'A_MISMATCHES': a_hit['mismatches'],
                                'A_BITSCORE': a_hit['bit_score'],
                                'B_LENGTH': get_gene_length('B', a_hit['subject']),
                                'B_COORDS': f"{b_hit['query_start']}-{b_hit['query_end']}",
                                'B_EVAL': b_hit['eval'],
                                'B_ALIGNMENT_LENGTH': b_hit['alignment_length'],
                                'B_MISMATCHES': b_hit['mismatches'],
                                'B_BITSCORE': b_hit['bit_score']
                            }
                            ),
                            ignore_index=True)

                # else:
                #     print(f'{colors.bcolors.RED}Not besties{colors.bcolors.END}')

    for locus_tag, a_hit in best_dict['A'].items():
        # print(f"---{locus_tag}-==")
        if a_hit['already_paired'] == False:
            #print(f'{colors.bcolors.RED}No reciprocal hit{a_hit}{colors.bcolors.END}')
            df = df.append(
                pd.Series(data={
                    'RECIPROCAL': "A",
                    'A_LOCUS': locus_tag,
                    'B_LOCUS': a_hit['subject'],
                    'A_LENGTH': get_gene_length('A', locus_tag),
                    'A_COORDS': f"{a_hit['query_start']}-{a_hit['query_end']}",
                    'A_EVAL': a_hit['eval'],
                    'A_ALIGNMENT_LENGTH': a_hit['alignment_length'],
                    'A_MISMATCHES': a_hit['mismatches'],
                    'A_BITSCORE': a_hit['bit_score'],
                    'B_LENGTH': get_gene_length('B', a_hit['subject'])
                }
                ),
                ignore_index=True)

    for locus_tag, b_hit in best_dict['B'].items():
        # print(f"---{locus_tag}-==")
        if b_hit['already_paired'] == False:
            #print(f'{colors.bcolors.RED}No reciprocal hit{b_hit}{colors.bcolors.END}')
            df = df.append(
                pd.Series(data={
                    'RECIPROCAL': "B",
                    'A_LOCUS': b_hit['subject'],
                    'B_LOCUS': locus_tag,
                    'A_LENGTH': get_gene_length('A', b_hit['subject']),
                    'B_LENGTH': get_gene_length('B', locus_tag),
                    'B_COORDS': f"{b_hit['query_start']}-{b_hit['query_end']}",
                    'B_EVAL': b_hit['eval'],
                    'B_ALIGNMENT_LENGTH': b_hit['alignment_length'],
                    'B_MISMATCHES': b_hit['mismatches'],
                    'B_BITSCORE': b_hit['bit_score']
                }
                ),
                ignore_index=True)

    # now find all genes that didn't even have a blast hit
    for locus_tag, gene in gene_dict['A'].items():
        if locus_tag not in best_dict['A']:
            df = df.append(
                pd.Series(data={
                    'RECIPROCAL': "A",
                    'A_LOCUS': locus_tag,
                    'B_LOCUS': "None",
                    'A_LENGTH': get_gene_length('A', locus_tag),

                }
                ),
                ignore_index=True)

    for locus_tag, gene in gene_dict['B'].items():
        if locus_tag not in best_dict['B']:
            df = df.append(
                pd.Series(data={
                    'RECIPROCAL': "B",
                    'A_LOCUS': "None",
                    'B_LOCUS': locus_tag,
                    'B_LENGTH': get_gene_length('B', locus_tag),

                }
                ),
                ignore_index=True)

    df.to_csv(args.out, sep="\t", index=False)
