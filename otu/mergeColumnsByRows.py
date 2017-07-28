# given two tables, it will merge the columns of the two tables into the same rows by row name

import sys

merged = {}

f1 = open(sys.argv[1],'r')
file1Lines = f1.readlines()
f1.close()

f2 = open(sys.argv[2],'r')
file2Lines = f2.readlines()
f2.close()

# get headers
headerBoth = []
header1 = file1Lines.pop(0).rstrip()
split = header1.split("\t")
header1 = split
for columnName in split:
    headerBoth.append(columnName)

header2 = file2Lines.pop(0).rstrip()
split = header2.split("\t")
header2 = split
for columnName in split:
    headerBoth.append(columnName)

# file 1

for line in file1Lines:
    line = line.rstrip()
    lineSplit = line.split("\t")
    rowName = lineSplit.pop(0)

    if rowName not in merged:

        merged[rowName] = {}
        for columnName in headerBoth:
            merged[rowName][columnName] = 0

    counter = 0
    while counter < len(header1):
        merged[rowName][header1[counter]] = lineSplit[counter]
        counter += 1

# file 2

for line in file2Lines:
    line = line.rstrip()
    lineSplit = line.split("\t")
    rowName = lineSplit.pop(0)

    if rowName not in merged:

        merged[rowName] = {}
        for columnName in headerBoth:
            merged[rowName][columnName] = 0

    counter = 0
    while counter < len(header2):
        merged[rowName][header2[counter]] = lineSplit[counter]
        counter += 1

# print combined

print(*headerBoth, sep = "\t")

for rowName in merged:
    print(rowName, end = "")
    for columnName in headerBoth:
        print("\t", merged[rowName][columnName], end = "")
    print()
