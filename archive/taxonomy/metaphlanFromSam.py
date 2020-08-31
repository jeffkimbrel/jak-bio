import os
import sys

def systemCall(command):
    print("\nSYSTEM CALL: "+command)
    os.system(command)

fileName = sys.argv[1]
handle = fileName.split("/")[-1].split(".")[0]

systemCall('python ~/Dropbox/scripts/biobakery-metaphlan2-c3fb65390c21/metaphlan2.py '+fileName+' --input_type sam --nproc 4 > '+handle+'.txt')
systemCall('python ~/Dropbox/scripts/bio/taxonomy/metaphlan2_2_edgeList.py -k BAVE -c 0.001 -f '+handle+'.txt > '+handle+'.edgeList.txt')
