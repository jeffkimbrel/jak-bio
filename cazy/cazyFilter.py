import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description = 'filters hmmsearch output based on the requirements set forth in the dbcan hmmscan-parser.sh file')
# (http://csbl.bmb.uga.edu/dbCAN/download/readme.txt)

parser.add_argument('-f', '--file', help = "Output from hmmsearch", required = True)
parser.add_argument('-l', '--length', default = 80, help = "Minimum alignment length for a 'good' hit" )
parser.add_argument('-c', '--covered', default = 0.3, help = "Minimum fraction of hmm covered for a 'bad' hit" )
parser.add_argument('--subfamily', '-s', action='store_true', help = 'Keep subfamilies' )
parser.add_argument('--verbose', '-v', action='store_true', help='Print original line' )

args = parser.parse_args()
args.length = int(args.length)
args.covered = float(args.covered)

lines = [line.strip() for line in open(args.file)]

for line in lines:
    if not line.lstrip().startswith('#'):

        passed = False

        split = line.split()
        model = split[3]

        model = model.split(".")[0]

        if args.subfamily == False:
            model = model.split("_")[0]

        evalue = float(split[6])
        alignmentLength = int(split[18]) - int(split[17]) + 1
        hmmCovered = (int(split[16]) - int(split[15]) + 1) / int(split[5])

        if alignmentLength >= args.length:
            if evalue <= 1e-05:
                passed = True
        elif evalue <= 1e-03:
            if hmmCovered >= args.covered:
                passed = True

        if passed == True:
            if args.verbose == True:
                print(line)
            else:
                print(split[0], model, sep = "\t")
