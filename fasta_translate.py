import sys
import argparse

from Bio import SeqIO
from Bio.Data import CodonTable

from dnachisel.biotools import reverse_translate
from collections import Counter

from jakomics import colors
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-f', '--fasta',
                    required=True)

parser.add_argument('--frames',
                    nargs="*",
                    default=[1])

parser.add_argument('--reverse',
                    action='store_true',
                    help='Reverse translate a fasta file of proteins')

args = parser.parse_args()

###

def report_problem_aa(seq_record):
    standard_aa = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    problems = False

    aa_freq = Counter(str(seq_record.seq))
    aa_prob_freq = {k:aa_freq[k] for k in aa_freq.keys() if k not in standard_aa}

    if len(aa_prob_freq) > 0:
        print(f"{colors.bcolors.YELLOW}WARNING\tnon-standard AAs\t{aa_prob_freq}\t{seq_record.description}{colors.bcolors.END}", file = sys.stderr)
        problems = True

    if problems:
        return(1)
    else:
        return(0)

def fix_ambiguous_aa(seq, add_stop):

    '''
    The pyrrolysine and selenocysteine can be changed at the codon table step in the main code,
    but other IUPAC ambiguous characters need to just be changed to one of the 
    '''

    replacements = {
        "X": "O", 
        "B": 'D', # Aspartic acid or Asparagine
        "Z": 'E'  # Glutamic acid or Glutamine
    }

    seq = str(seq).upper() + "*"

    for aa, codon in replacements.items():
        seq = seq.replace(aa, codon)

    return(seq)


###

if __name__ == "__main__":
    jak_utils.header()

    if args.reverse == False:
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
    
    else:
        CodonTable.unambiguous_dna_by_name['Standard'].back_table['U'] = 'TGA'
        CodonTable.unambiguous_dna_by_name['Standard'].back_table['O'] = 'TAG'
        CodonTable.unambiguous_dna_by_name['Standard'].forward_table['TGA'] = 'U'
        CodonTable.unambiguous_dna_by_name['Standard'].forward_table['TAG'] = 'O'

        failed = 0
        warnings = 0

        for seq_record in SeqIO.parse(args.fasta, "fasta"):

            # report problems
            warnings = warnings + report_problem_aa(seq_record)

            # fix problematic characters
            fixed_seq = fix_ambiguous_aa(seq_record.seq, add_stop = True)

            try:
                s = reverse_translate(fixed_seq, table="Standard")
                print(f">{seq_record.description}")
                print(s)
            except KeyError as e:
                print(f"{colors.bcolors.RED}ERROR - {e} not found: {seq_record.id}{colors.bcolors.END}", file=sys.stderr)
                failed += 1
                continue

print(f"{colors.bcolors.YELLOW}Warnings: {warnings}{colors.bcolors.END}", file=sys.stderr)
print(f"{colors.bcolors.RED}Failures: {failed}{colors.bcolors.END}", file=sys.stderr)