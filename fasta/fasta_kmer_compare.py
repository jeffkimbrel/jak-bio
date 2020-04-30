from itertools import combinations
import os
import argparse
import subprocess

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='')

parser.add_argument('-d', '--directory',
                    required=True)

parser.add_argument('-k', '--kmer',
                    default=[21],
                    nargs='*',
                    help='kmer size')

args = parser.parse_args()

# FUNCTIONS ###################################################################


def get_genome_pairs():
    fileList = []

    for fileName in os.listdir(args.directory):
        if fileName.endswith('fna'):
            fileList.append(fileName)
    return list(combinations(fileList, 2))


def kat_comp(pair, k):
    g1_name = os.path.splitext(os.path.basename(pair[0]))[0]
    g2_name = os.path.splitext(os.path.basename(pair[1]))[0]

    output_file = 'kat/' + g1_name + '_' + g2_name + '_' + str(k)

    command = 'kat comp -t 8 -m ' + str(k) + ' -o ' + \
        output_file + ' ' + pair[0] + ' ' + pair[1]
    # os.system(command)

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    stats = stdout.decode().split("\n")

    # total kmers
    g1_t = stats[24].split(": ")[1]
    g2_t = stats[25].split(": ")[1]

    # distinct kmers
    g1_d = stats[28].split(": ")[1]
    g2_d = stats[29].split(": ")[1]

    # distinct kmers in genome only
    g1_do = stats[36].split(": ")[1]
    g2_do = stats[37].split(": ")[1]

    # shared kmers
    shared = stats[42].split(": ")[1]

    print(k, g1_name, g2_name, g1_t, g2_t, g1_d, g2_d, g1_do, g2_do, shared, sep="\t")


# LOOP ########################################################################
genome_pairs = get_genome_pairs()

for pair in genome_pairs:
    for k in args.kmer:
        kat_comp(pair, k)
