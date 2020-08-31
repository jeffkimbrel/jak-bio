import sys
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# give it a fasta file and a number (representing a barcode length) and it will find and count all instances at the 5' end

barcodeSize = sys.argv[1]
barcodes = {}

for record in SeqIO.parse(sys.argv[1], "fasta"):
    barcode = record.seq[:int(sys.argv[2])]
    barcodes[barcode] = barcodes.get(barcode, 0) + 1

for barcode in barcodes:
    print(barcode, barcodes[barcode], sep = "\t")
