import argparse
import os
from Bio import SeqIO

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = 'X')
parser.add_argument('-i1', '--index1',  help = "Index 1 fastq",     required = True)
parser.add_argument('-i2', '--index2',  help = "Index 2 fastq",     required = True)
parser.add_argument('-r1', '--read1',   help = "Read 1 fastq",      required = True)
parser.add_argument('-r2', '--read2',   help = "Read 2 fastq",      required = True)
parser.add_argument('-m',  '--mapping', help = "Mapping File",      required = True)
parser.add_argument('-b',  '--buffer',  help = "read buffer size",  required = False, default = 10000)
parser.add_argument('-f',  '--folder',  help = "Output Folder",     required = False, default = "demux")

args = parser.parse_args()

if not os.path.exists(args.folder):
    os.makedirs(args.folder)

total = 0

## CLASS #######################################################################

samples = {}
class Sample:
    def __init__(self, sample, barcode):
        self.sample = sample
        self.barcode = barcode
        self.file1 = args.folder + "/" + sample + "_R1.fastq"
        self.file2 = args.folder + "/" + sample + "_R2.fastq"
        self.read1 = []
        self.read2 = []
        self.count = 0

        # create empty files
        open(self.file1, 'w').close()
        open(self.file2, 'w').close()

    def addPair(self, read1, read2):
        self.read1.append(read1)
        self.read2.append(read2)
        self.count += 1

    def write_to_file(self, final):
        if final == True or len(self.read1) > args.buffer:
                output_handle = open(self.file1, "a")
                SeqIO.write(self.read1, output_handle, "fastq")
                self.read1 = []

                output_handle = open(self.file2, "a")
                SeqIO.write(self.read2, output_handle, "fastq")
                self.read2 = []

## FUNCTIONS ###################################################################

def createUndetermined():
    samples['Undetermined'] = Sample('Undetermined', 'Undetermined')

def addSamples():
    file = [line.strip() for line in open(args.mapping)]
    for line in file:
        if not line.startswith('#'):
            split = line.split("\t")
            samples[split[1]] = Sample(split[0], split[1])

def readFiles():

    total = 0
    errors = 0

    for read1, read2, index1, index2 in zip(SeqIO.parse(args.read1,  "fastq"),
                                            SeqIO.parse(args.read2,  "fastq"),
                                            SeqIO.parse(args.index1, "fastq"),
                                            SeqIO.parse(args.index2, "fastq")):

        headers_agree = False

        if len(set({read1.id, read2.id, index1.id, index2.id})) == 1:
            headers_agree = True
        else:
            errors += 1

        barcode = index1.seq + index2.seq

        if not barcode in samples:
            barcode = 'Undetermined'

        if headers_agree == True:

            samples[barcode].addPair(read1, read2)
            samples[barcode].write_to_file(False)

        total += 1
        if total % args.buffer == 0:
            print("Processed " + str(total) + " total pairs with " + str(errors) + " pairing errors            \r", end = '')

    return(total, errors)

def main():
    createUndetermined()
    addSamples()

    total, errors = readFiles()

    print("Processed " + str(total) + " total pairs with " + str(errors) + " pairing errors")

    # a final write
    for barcode in samples:
        samples[barcode].write_to_file(True)

## RUN #########################################################################

if __name__== "__main__":
    main()
