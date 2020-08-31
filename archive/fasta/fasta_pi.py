from Bio import SeqIO
from Bio.SeqUtils import ProtParam
#import os.path
import argparse

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='Give it a multifasta file of amino acids, and it will draw a histogram of the theoretical isoelectric points (pIs). This version uses R and ggplot to produce the plot.')

parser.add_argument('-f', '--file',
	help='Fasta file with aa sequences')
parser.add_argument('-v', '--version',
	action='version',
	version='%(prog)s 2.0')
args = parser.parse_args()

## MAIN ########################################################################

for seq_record in SeqIO.parse(args.file, "fasta"):
	PP = ProtParam.ProteinAnalysis(str(seq_record.seq))
	print(seq_record.id, PP.isoelectric_point(), sep = "\t")

	#pI_list.append(PP.isoelectric_point())
