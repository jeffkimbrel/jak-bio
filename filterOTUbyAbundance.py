import sys

minAbundance = 2
minSample = 2

f = open(sys.argv[1],'r')
print(f.readlines()[0])
f.close()

f = open(sys.argv[1],'r')
lines = f.readlines()[1:]
f.close()

for line in lines:
    line = line.strip()
    split = line.split("\t")
    passed = 0

    itersplit = iter(split)
    next(itersplit)

    for value in itersplit:
        if int(value) >= minAbundance:
            passed = passed + 1

    if passed >= minSample:
        print(line)
