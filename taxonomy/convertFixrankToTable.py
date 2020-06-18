import sys
import argparse

# Arguments
parser = argparse.ArgumentParser(description='XXX')

parser.add_argument('-f', '--file',
                    help="RDP fixrank file",
                    required=True)

parser.add_argument('-t', '--threshold',
                    help="RDP confidence threshold",
                    type=float,
                    default=0.8,
                    required=True)

args = parser.parse_args()

if args.threshold > 1:
    args.threshold /= 100

###


def p2f(x):
    return float(x.strip('%'))/100


f = open(args.file, 'r')
lines = f.readlines()[7:]
f.close()

# header
print("", "domain", "phylum", "class", "order", "family", "genus", sep="\t")


for line in lines:
    line = line.strip()
    split = line.split(";")

    # print(line)

    print(split[0], end="\t")

    # last success
    lastSuccess = "root"

    # domain
    domain = ""
    if p2f(split[3]) >= args.threshold:
        domain = split[2]
        domain = domain.replace('"', '')
        lastSuccess = domain
    else:
        domain = "unclassified_" + lastSuccess

    print(domain, end="\t")

    # phylum
    phylum = ""
    if p2f(split[5]) >= args.threshold:
        phylum = split[4]
        phylum = phylum.replace('"', '')
        lastSuccess = phylum
    else:
        phylum = "unclassified_" + lastSuccess

    print(phylum, end="\t")

    # class
    classC = ""
    if p2f(split[7]) >= args.threshold:
        classC = split[6]
        classC = classC.replace('"', '')
        lastSuccess = classC
    else:
        classC = "unclassified_" + lastSuccess

    print(classC, end="\t")

    # order
    order = ""
    if p2f(split[9]) >= args.threshold:
        order = split[8]
        order = order.replace('"', '')
        lastSuccess = order
    else:
        order = "unclassified_" + lastSuccess

    print(order, end="\t")

    # family
    family = ""
    if p2f(split[11]) >= args.threshold:
        family = split[10]
        family = family.replace('"', '')
        lastSuccess = family
    else:
        family = "unclassified_" + lastSuccess

    print(family, end="\t")

    # genus
    genus = ""
    if p2f(split[13]) >= args.threshold:
        genus = split[12]
        genus = genus.replace('"', '')
        lastSuccess = genus
    else:
        genus = "unclassified_" + lastSuccess

    print(genus)

    # print()
