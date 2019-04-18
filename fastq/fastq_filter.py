import argparse
import os
import sys
import subprocess

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(description = 'Takes all fastq pairs in -d and quality filters and kmer checks, writes to -o. Does not trim.')

parser.add_argument('-d', '--directory',
    help = "Directory with fastq files",
    required = True)
parser.add_argument('-o', '--out',
    default = "fastq_filter/",
    help = "Output Folder")

parser.add_argument('--quality', '-q',
    action = 'store_true',
    help = 'quality filter' )

parser.add_argument('--contaminants', '-c',
    action = 'store_true',
    help = 'filter contaminants' )

parser.add_argument('--test', '-t',
    action = 'store_true',
    help = 'Test commands without running them' )

args = parser.parse_args()

args.directory = os.path.abspath(args.directory)
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

## FUNCTIONS ###################################################################

def identify_pairs(path):

    # identify fastq pairs with the same sample name (text before first underscore)

    d = os.listdir(path)
    pairs = {}

    for fileName in d:
        if fileName.endswith('fastq') or fileName.endswith('fastq.gz'):
            split = fileName.split("_")
            sample = split[0]

            if not sample in pairs:
                pairs[sample] = {"R1" : [], "R2" : []}

            if "R1" in split:
                pairs[sample]['R1'].append(fileName)

            if "R2" in split:
                pairs[sample]['R2'].append(fileName)

    return(pairs)

def verify_pairs(pairs):

    # makes sure each pair has exactly on set of forward and reverse reads.

    for sample in pairs:
        if len(pairs[sample]["R1"]) != 1 or len(pairs[sample]["R2"]) != 1:
            print(sample + " has a read pair error!! Stopping.")
            sys.exit()

    print("\n*** All fastq file names in " + args.directory + " appear to be correct! Nice!!\n")

def extract_stats(lines):
    stats = {"input" : 0, "low" : 0, "contamination" : 0, "removed" : 0, "remain" : 0}

    for line in lines:
        if line.startswith("Input:"):
            stats["input"] = line.split()[1]
        elif line.startswith("Contaminants:"):
            stats["contamination"] = line.split()[1]
        elif line.startswith("Low quality discards:"):
            stats["low"] = line.split()[3]
        elif line.startswith("Total Removed:"):
            stats["removed"] = line.split()[2]
        elif line.startswith("Result:"):
            stats["remain"] = line.split()[1]

    return(stats)

def filter_pair(sample, pair):
    call = 'bbduk.sh in1=' + \
    args.directory + "/" + pair['R1'][0] + \
    ' in2=' + \
    args.directory + "/" + pair['R2'][0] + \
    ' out1=' + \
    args.out + "/" + sample + "_R1.filtered.fastq.gz" + \
    ' out2=' + \
    args.out + "/" + sample + "_R2.filtered.fastq.gz"

    if args.quality == True:
        call += " maq=20"

    if args.contaminants == True:
        call += " ref=" + os.path.dirname(os.path.abspath(sys.argv[0])) + "/contam_seqs.fa k=31 hdist=1 "

    print("---\n$> " + call)

    if args.test == True:
        p1 = subprocess.Popen(call, shell = True, stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = p1.communicate()
        err = err.decode()
        lines = err.split('\n')
        stats = extract_stats(lines)
        return(stats)
    else:
        return("TEST")

def format_results(results):
    print("SAMPLE", "READS", "CONTAMINANTS", "LOW-Q", "REMOVED", "REMAIN", sep = " | ")
    print("---", "---", "---", "---", "---", "---", sep = " | ")

    for sample in results:
        print(sample, results[sample]["input"], results[sample]["contamination"], results[sample]["low"], results[sample]["removed"], results[sample]["remain"], sep = " | ")


## MAIN ########################################################################

pairs = identify_pairs(args.directory)

verify_pairs(pairs)

results = {}

for sample in sorted(pairs.keys()):
    results[sample] = filter_pair(sample, pairs[sample])

if args.test == False:
    format_results(results)
