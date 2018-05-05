import sys

cutoff = 0.8 # 0.5 (50%) is recommended for reads less than 300bp. However, this might have changed to 250bp at some recent version of RDP: https://rdp.cme.msu.edu/classifier/class_help.jsp#conf

def p2f(x):
    return float(x.strip('%'))/100

f = open(sys.argv[1],'r')
lines = f.readlines()[7:]
f.close()

# header
print("", "domain", "phylum", "class", "order", "family", "genus", sep = "\t")


for line in lines:
    line = line.strip()
    split = line.split(";")

    #print(line)

    print(split[0],end = "\t")

    # last success
    lastSuccess = "root"

    # domain
    domain = ""
    if p2f(split[3]) >= cutoff:
        domain = split[2]
        domain = domain.replace('"', '')
        lastSuccess = domain
    else:
        domain = "unclassified_" + lastSuccess

    print(domain, end = "\t")

    # phylum
    phylum = ""
    if p2f(split[5]) >= cutoff:
        phylum = split[4]
        phylum = phylum.replace('"', '')
        lastSuccess = phylum
    else:
        phylum = "unclassified_" + lastSuccess

    print(phylum, end = "\t")

    # class
    classC = ""
    if p2f(split[7]) >= cutoff:
        classC = split[6]
        classC = classC.replace('"', '')
        lastSuccess = classC
    else:
        classC = "unclassified_" + lastSuccess

    print(classC, end = "\t")


    # order
    order = ""
    if p2f(split[9]) >= cutoff:
        order = split[8]
        order = order.replace('"', '')
        lastSuccess = order
    else:
        order = "unclassified_" + lastSuccess

    print(order, end = "\t")

    # family
    family = ""
    if p2f(split[11]) >= cutoff:
        family = split[10]
        family = family.replace('"', '')
        lastSuccess = family
    else:
        family = "unclassified_" + lastSuccess

    print(family, end = "\t")

    # genus
    genus = ""
    if p2f(split[13]) >= cutoff:
        genus = split[12]
        genus = genus.replace('"', '')
        lastSuccess = genus
    else:
        genus = "unclassified_" + lastSuccess

    print(genus, end = "\t")

    print()
