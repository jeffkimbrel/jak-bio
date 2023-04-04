#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json

from jakomics import colors, utilities
import jak_utils
jak_utils.header()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Filters fasta files, writing to a new file appended with ".ff"')

parser.add_argument('-f', '--file',
                    help="Paths to individual fasta files",
                    required=True)

parser.add_argument('-m', '--mode',
                    help="id or description mode",
                    required=False,
                    default="description")


parser.add_argument('--remove_identical_sequences',
                    action='store_true',
                    help='Remove stop codon asterisks')

args = parser.parse_args()

if args.mode not in ["description", "id"]:
    print(f"{colors.bcolors.RED}ERROR: --mode (-m) must be either 'description' or 'id' {colors.bcolors.END}", file=sys.stderr)
    sys.exit()

# FUNCTIONS

def get_header(seq):
    # chosen based off the --mode selection

    if args.mode == 'description':
        header = seq.description
    elif args.mode == "id":
        header = seq.id
    return(header)

def validate_records(seqs): 

    deduped_records = {}

    for seq in seqs:

        header = get_header(seq)

        if header not in deduped_records:
            deduped_records[header] = [str(seq.seq)]
        else:
            deduped_records[header].append(str(seq.seq))

    good_records = 0
    bad_records = {'same_sequence': {},
                   'diff_sequence': {}}


    for header, records in deduped_records.items():
        if len(records) == 1:
            good_records += 1

        if len(records) > 1:
            
            # check if the sequences are identical
            if all(element == records[0] for element in records):
                bad_records['same_sequence'][header] = records
            else:
                bad_records['diff_sequence'][header] = records

    print(f"There are {len(seqs)} total records")
    print(f"Of these, {good_records} look OK")
    print(f"{len(bad_records['same_sequence'])} records are found more than once, but each record has the same sequence")
    print(f"{len(bad_records['diff_sequence'])} records are found more than once, and there is more than one sequence")

    with open("out.json", 'w') as outfile:
        json.dump(bad_records, outfile, indent=2)






# MAIN ########################################################################

#pbar = tqdm(total=len(fasta_files), desc="Finished", unit=" fasta files")

if __name__ == "__main__":

    seqs = []
    for seq_record in SeqIO.parse(args.file, "fasta"):
        seqs.append(seq_record)

    validate_records(seqs)
