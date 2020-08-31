import sys

level = sys.argv[3] # taxonomic level of comparison

f = open(sys.argv[1],'r')
linesOld = f.readlines()[1:]
f.close()

f = open(sys.argv[2],'r')
linesNew = f.readlines()[1:]
f.close()

# process old
dictOld = {}

for line in linesOld:
    line = line.strip()
    split = line.split("\t")

    dictOld[split[0]] = {"domain" : split[1].replace('"', ''),
                         "phylum" : split[2].replace('"', ''),
                         "class"  : split[3].replace('"', ''),
                         "order"  : split[4].replace('"', ''),
                         "family" : split[5].replace('"', ''),
                         "genus"  : split[6].replace('"', '')
                         }

# process new
dictNew = {}

for line in linesNew:
    line = line.strip()
    split = line.split("\t")

    dictNew[split[0]] = {"domain" : split[1].replace('"', ''),
                         "phylum" : split[2].replace('"', ''),
                         "class"  : split[3].replace('"', ''),
                         "order"  : split[4].replace('"', ''),
                         "family" : split[5].replace('"', ''),
                         "genus"  : split[6].replace('"', '')
                         }

# compare

changes = {}
agreement = {"yes" : 0, "no" : 0}
unclassified = {"old" : 0, "new" : 0}


for RSV in dictOld:
    old = dictOld[RSV][level]
    if "unclassified" in old:
        unclassified["old"] += 1

    new = dictNew[RSV][level]
    if "unclassified" in new:
        unclassified["new"] += 1

    if old != new:
        agreement["no"] += 1
        if old in changes:
            changes[old].append(new)
        else:
            changes[old] = [new]
    else:
        agreement["yes"] += 1

for change in changes:
    changes[change] = list(set(changes[change]))
    print(change, ": ", sep = "", end = "")
    print(*changes[change], sep = ", ")

print(agreement)
print(unclassified)
