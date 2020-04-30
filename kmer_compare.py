from Bio import SeqIO
from Bio.Seq import Seq
import math

kmer_length = 31
window_size = 250

# FUNCTIONS ###################################################################


def get_ref_kmers():
    reference_seq = ""
    reference_kmers = {}

    for seq_record in SeqIO.parse("/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/CSP1939/isolates/Algoriphagus_sp_ARW1R1/JGI/IMG/2747842515.fna", "fasta"):
        reference_seq = reference_seq + "|" + str(seq_record.seq)

    for i in range(len(reference_seq)):
        kmer = reference_seq[i:i+kmer_length]

        if "|" not in kmer:
            rc = str(Seq(kmer).reverse_complement())

            if kmer in reference_kmers:
                reference_kmers[kmer]['ref'].append(i)
            elif rc in reference_kmers:
                reference_kmers[rc]['ref'].append(i)
            else:
                reference_kmers[kmer] = {'ref': [i], "comp": []}

    return reference_kmers


def get_comp_kmers(reference_kmers):
    compare_seq = ""

    for seq_record in SeqIO.parse("/Users/kimbrel1/Dropbox/LLNL/Projects/Biofuels_SFA/CSP1939/isolates/Marinobacter_sp_PT3-2/JGI/IMG/Ga0264245_total_v1.fna", "fasta"):
        compare_seq = compare_seq + "|" + str(seq_record.seq)

    for i in range(len(compare_seq)):
        kmer = compare_seq[i:i+kmer_length]

        if "|" not in kmer:
            rc = str(Seq(kmer).reverse_complement())

            if kmer in reference_kmers:
                reference_kmers[kmer]['comp'].append(i)
            elif rc in reference_kmers:
                reference_kmers[rc]['comp'].append(i)
            else:
                # kmers found in compare, but not in reference (should be rare?)
                reference_kmers[kmer] = {'ref': [], "comp": [i]}

    return reference_kmers


reference_kmers = get_ref_kmers()
reference_kmers = get_comp_kmers(reference_kmers)

reference_sums = {}
compare_sums = {}

for kmer in reference_kmers:

    for pos in reference_kmers[kmer]['ref']:
        if len(reference_kmers[kmer]['ref']) > 0:
            window = math.floor(pos / window_size)

            reference_sums[window] = reference_sums.get(
                window, 0) + len(reference_kmers[kmer]['ref'])
            compare_sums[window] = compare_sums.get(
                window, 0) + len(reference_kmers[kmer]['comp'])


for window in reference_sums:
    print(window * window_size, reference_sums[window] /
          window_size, compare_sums[window] / window_size, sep="\t")
