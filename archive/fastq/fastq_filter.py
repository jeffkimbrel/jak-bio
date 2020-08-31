import argparse
import os
import sys
import subprocess

## OPTIONS #####################################################################

parser = argparse.ArgumentParser(
    description='Takes all fastq pairs in -d and quality filters and kmer checks, writes to -o. Does not trim.')

parser.add_argument('-d', '--directory',
                    help="Directory with fastq files",
                    required=True)

parser.add_argument('-o', '--out',
                    default="fastq_filter/",
                    help="Output Folder")

parser.add_argument('-q', '--quality',
                    default=None,
                    help="Minimum average quality")

parser.add_argument('--contaminants', '-c',
                    action='store_true',
                    help='filter contaminants')

args = parser.parse_args()

args.directory = os.path.abspath(args.directory)
args.out = os.path.abspath(args.out)

if not os.path.exists(args.out):
    os.makedirs(args.out)

## CLASSES #####################################################################

samples = {}


class Sample:
    def __init__(self, name):
        self.name = name
        self.R1 = []
        self.R2 = []
        self.valid_file_pairs = True
        self.ordered = False

    def add_file(self, file):
        if "R1" in file:
            self.R1.append(args.directory + "/" + file)

        if "R2" in file:
            self.R2.append(args.directory + "/" + file)

    def validate_file_pairs(self):
        if len(self.R1) != 1 or len(self.R2) != 1:
            self.valid_file_pairs = False

        for R1 in self.R1:
            print("R1_in:", R1, sep="\t")
        for R2 in self.R2:
            print("R2_in:", R2, sep="\t")

        if self.valid_file_pairs == False:
            print("\n***  ERROR  ***\n" + self.name +
                  " doesn't have proper file names. Each sample should have exactly one R1 and one R2 file. \n*** EXITING ***\n")
            sys.exit()
        else:
            print("Proper file pairs have been found.")

        return(self.valid_file_pairs)

    def verify_read_pairs(self):
        call = 'reformat.sh in1=' + self.R1[0] + ' in2=' + self.R2[0] + ' verifypaired=t'
        lines = system_call(call)

        if "Names appear to be correctly paired." in lines:
            self.ordered = True
            print("Reads are in the same order in both files.")
        else:
            print("\n***  ERROR  ***\nThe read pairs are not in the same order in both files for " +
                  self.name + "\n*** EXITING ***\n")
            sys.exit()

    def set_out_names(self):
        self.R1_out = args.out + "/" + self.name + "_filtered.R1.fastq.gz"
        self.R2_out = args.out + "/" + self.name + "_filtered.R2.fastq.gz"

        print("R1_out:", self.R1_out, sep="\t")
        print("R2_out:", self.R2_out, sep="\t")

    def filter_pair(self):
        call = 'bbduk.sh in1=' + self.R1[0] + ' in2=' + \
            self.R2[0] + ' out1=' + self.R1_out + ' out2=' + self.R2_out

        if args.quality != None:
            call += " maq=" + str(args.quality)
            print("Quality filtering at maq=" + str(args.quality))
        if args.contaminants == True:
            call += " ref=" + \
                os.path.dirname(os.path.abspath(sys.argv[0])) + "/contam_seqs.fa k=31 hdist=1 "
            print("Removing contaminants from " +
                  os.path.dirname(os.path.abspath(sys.argv[0])) + "/contam_seqs.fa")

        # stats
        call += " stats=" + args.out + "/" + self.name + "_stats.txt "

        self.stats = extract_stats(system_call(call))

    def print_verbose_stats(self):

        print("\nReads in:", self.stats["input"], sep="\t\t")
        print("Order Verified:", self.ordered, sep="\t\t")
        print("Contaminants Removed:", self.stats["contamination"], cp(
            self.stats["contamination"], self.stats["input"]), sep="\t")
        print("Low-Quality Removed:", self.stats["lowQC"],
              cp(self.stats["lowQC"], self.stats["input"]), sep="\t")
        print("Total Removed:\t", self.stats["removed"], cp(
            self.stats["removed"], self.stats["input"]), sep="\t")
        print("Remaining Reads:", self.stats["remain"], cp(
            self.stats["remain"], self.stats["input"]), sep="\t")
        print()

## FUNCTIONS ###################################################################


def system_call(call):
    p1 = subprocess.Popen(call, shell=True, stdin=None,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p1.communicate()
    err = err.decode()
    lines = err.split('\n')
    return(lines)


def cp(n, d):  # convert percentage
    return(str("{0:.6g}".format(100 * int(n) / int(d))) + "%")


def add_files():
    files = os.listdir(args.directory)
    for file in sorted(files):
        if file.endswith('fastq') or file.endswith('fastq.gz'):
            split = file.split("_")
            sample = split[0]

            if sample not in samples:
                samples[sample] = Sample(sample)

            samples[sample].add_file(file)


def extract_stats(lines):
    stats = {"input": 0, "lowQC": 0, "contamination": 0, "removed": 0, "remain": 0}

    for line in lines:
        if line.startswith("Input:"):
            stats["input"] = line.split()[1]
        elif line.startswith("Contaminants:"):
            stats["contamination"] = line.split()[1]
        elif line.startswith("Low quality discards:"):
            stats["lowQC"] = line.split()[3]
        elif line.startswith("Total Removed:"):
            stats["removed"] = line.split()[2]
        elif line.startswith("Result:"):
            stats["remain"] = line.split()[1]

    return(stats)


def print_all_stats():
    print("SAMPLE", "READS", "ORDERED", "CONTAMINANTS", "LOWQ", "REMOVED", "REMAIN", sep="\t")
    for sample in samples:
        print(sample, samples[sample].stats["input"], samples[sample].ordered, samples[sample].stats["contamination"],
              samples[sample].stats["lowQC"], samples[sample].stats["removed"], samples[sample].stats["remain"], sep="\t")


def main():
    add_files()

    for number, sample in enumerate(samples):
        print("--- " + sample + " (" + str(number + 1) + "/" +
              str(len(samples)) + ") ----------------------------------")
        samples[sample].validate_file_pairs()
        samples[sample].verify_read_pairs()
        samples[sample].set_out_names()
        samples[sample].filter_pair()
        samples[sample].print_verbose_stats()

    print_all_stats()


main()
