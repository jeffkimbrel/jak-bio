import sys
import argparse

parser = argparse.ArgumentParser(description = 'Takes two datasets and merges by a common ID')

parser.add_argument('-f1', '--file1',
    help = "First file name",
    required = True)

parser.add_argument('-f2', '--file2',
    help = "Second file name",
    required = True)

parser.add_argument('-i1', '--id1',
    help = "First file id column (1-based)",
    type=int,
    default = 1)

parser.add_argument('-i2', '--id2',
    help = "Second file id column (1-based)",
    type=int,
    default = 1)

parser.add_argument('-d1', '--data1',
    help = "First file data column (1-based)",
    type=int,
    default = 2)

parser.add_argument('-d2', '--data2',
    help = "Second file data column (1-based)",
    type=int,
    default = 2)

args = parser.parse_args()

args.id1 -= 1
args.id2 -= 1
args.data1 -= 1
args.data2 -= 1

## LOOP ########################################################################

merged = {}

## FILE 1

file1 = [line.strip() for line in open(args.file1)]

for line in file1:
    split = line.rstrip().split("\t")
    id = split[args.id1]
    data = split[args.data1]

    if id in merged:
        merged[id]['file1'] = data
    else:
        merged[id] = {'file1' : data}


## FILE 2

file2 = [line.strip() for line in open(args.file2)]

for line in file2:
    split = line.rstrip().split("\t")
    id = split[args.id2]
    data = split[args.data2]

    if id in merged:
        merged[id]['file2'] = data
    else:
        merged[id] = {'file2' : data}

for id in merged:
    print(id, merged[id])
