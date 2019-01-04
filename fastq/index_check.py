from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

import sys

# read sample data
barcodes = {}
barcodes['undetermined'] = {'sample' : 'undetermined', 'count' : 0}
file = [line.strip() for line in open('SFA-1&2_BF-OvG-SEQ_G_qME_time_series_reseqs_16S_mapping_file_corrected.txt')]
for line in file:
    if not line.startswith('#'):
        split = line.split("\t")
        barcodes[split[1]] = {'sample' : split[0], 'count' : 0}

sequences = {}

for rec in SeqIO.parse("test_i1.fastq", "fastq"):
#for rec in SeqIO.parse("Undetermined_S0_L001_I1_001.fastq", "fastq"):
    if rec.id not in sequences:
        sequences[rec.id] = {'forward' : rec.seq}
    else:
        print("ERROR", file = sys.stderr)

for rec in SeqIO.parse("test_i2.fastq", "fastq"):
#for rec in SeqIO.parse("Undetermined_S0_L001_I2_001.fastq", "fastq"):
    if rec.id not in sequences:
        print("ERROR", file = sys.stderr)
    else:
        sequences[rec.id]['reverse'] = rec.seq

for seq in sequences:
    concat = sequences[seq]['forward'] + sequences[seq]['reverse']

    if concat in barcodes:
        barcodes[concat]['count'] += 1
    else:
        barcodes['undetermined']['count'] += 1

for barcode in sorted(barcodes.keys()):
    print(barcode, barcodes[barcode]['sample'], barcodes[barcode]['count'], sep = "\t")
