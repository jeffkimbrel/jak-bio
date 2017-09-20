import os
import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'XXX')

parser.add_argument('-f', '--file', help = "Amino acid fasta file", required = True)
parser.add_argument('-n', '--name', help = "Run name", required = True)
parser.add_argument('-o', '--out', help = "Output folder", required = True)

args = parser.parse_args()

unsortedName = args.out + "/" + args.name + ".unsorted.txt"
sortedName = args.out + "/" + args.name + ".sorted.txt"
resultsName = args.out + "/" + args.name + ".dbcan.txt"


def systemCall(command):
    print("\n### SYSTEM CALL: "+command)
    os.system(command)

systemCall("pwd")

# run hmmsearch
systemCall("hmmsearch -o " + args.out + "/dbcan.log --domT 10 --domtblout "  + unsortedName + " ~/Dropbox/scripts/hmm/dbCAN-fam-HMMs.txt.v5.txt " + args.file )
systemCall("sort " + unsortedName + " > " + sortedName)
systemCall("rm " + unsortedName)
systemCall("python ~/Dropbox/scripts/bio/cazy/cazy_filter.py -f " + sortedName + " > " + resultsName)

systemCall("python ~/Dropbox/scripts/bio/cazy/cazy_summary.py -m -f " + resultsName + " -t b | column -t")
systemCall("python ~/Dropbox/scripts/bio/cazy/cazy_summary.py -m -f " + resultsName + " -t c | column -t")
#systemCall("python ~/Dropbox/scripts/bio/cazy/cazy_summary.py -m -f " + resultsName + " -t f | column -t")
