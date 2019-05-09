import gzip
from Bio import SeqIO
import argparse
import subprocess

from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser(description = '')

parser.add_argument('-f', '--files',
    help = "fastq files",
    nargs = '*',
    required = True)

args = parser.parse_args()

##

def get_md5(file):
    call = "md5 " + file
    p1 = subprocess.Popen(call, shell = True, stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = p1.communicate()
    out = out.decode()
    junk, md5 = out.split(" = ")
    return(md5.strip())

##

print("FILE", "MD5", "INSTRUMENT", "RUN", "FLOWCELL", "LANE", "READS", sep = " | ")
print("---", "---", "---", "---", "---", "---", "---", sep = " | ")

for file in args.files:

    info = {}

    md5 = get_md5(file)

    with gzip.open(file, "rt") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
        #for record in SeqIO.parse(handle, "fastq"):
            #instrument, run, flowcell, lane, tile, x, y = record.id.split(":")
            split = title.split(":")
            merge = split[0] + ":" + split[1] + ":" + split[2] + ":" + split[3] + ":" + md5

            if merge in info:
                info[merge] += 1
            else:
                info[merge] = 1

    for record in info:
        instrument, run, flowcell, lane, md5 = record.split(":")
        print(file, md5, instrument, run, flowcell, lane, info[record], sep = " | ")
