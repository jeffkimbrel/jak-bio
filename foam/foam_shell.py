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
bhName = args.out + "/" + args.name + ".bh.txt"
koName = args.out + "/" + args.name + ".ko.txt"
tempA = args.out + "/A.temp"
tempB = args.out + "/B.temp"
tempC = args.out + "/C.temp"

def systemCall(command):
    print("\n### SYSTEM CALL: "+command)
    os.system(command)

systemCall("pwd")

# run hmmsearch
systemCall("hmmsearch -o " + args.out + "/dbcan.log --domT 25 --cpu 4 --domtblout "  + unsortedName + " ~/Dropbox/scripts/FOAM/FOAM-hmm_rel1a.hmm " + args.file )
systemCall("sort " + unsortedName + " > " + sortedName)
systemCall("rm " + unsortedName)

systemCall("python ~/Dropbox/scripts/bio/foam/bmn-HMMerBestHit.py " + sortedName + " > " + tempA)
systemCall("python ~/Dropbox/scripts/bio/foam/cleanupBH.py " + tempA + " > " + bhName)

systemCall("awk '{print $4}' " + bhName + " > " + tempB)
systemCall("python ~/Dropbox/scripts/bio/foam/bmn-CountEachElement.py " + tempB + " > " + tempC)
systemCall("python ~/Dropbox/scripts/bio/foam/bmn-KOoneCount.py " + tempC + " | sed s/KO://g | sort -k 1 > " + koName)
