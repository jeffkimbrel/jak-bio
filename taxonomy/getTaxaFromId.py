import sys
import time
from Bio import Entrez
Entrez.email = 'kimbrel1@llnl.gov'

dictionary = {}

## Read in list
lines = [line.strip() for line in open(sys.argv[1])]

for line in list(lines):
    split = line.split("|")

    inputID = split[3]

    print(inputID,file=sys.stderr)

    dictionary[inputID] = {}

    #
    #time.sleep(1) # add some delay for NCBI sake
    handle = Entrez.esearch(db="protein", term=inputID, retmode="xml")
    record = Entrez.read(handle)
    accession = record["IdList"][0]

    #
    #time.sleep(1) # add some delay for NCBI sake
    handle = Entrez.esummary(db="protein", id=accession, retmode="xml")
    record = Entrez.read(handle)

    tax = record[0]['TaxId']

    for item in record[0]:
        dictionary[inputID][item] = record[0][item]

    #
    #time.sleep(1) # add some delay for NCBI sake
    handle = Entrez.efetch(db="taxonomy", id=str(tax), retmode="xml")
    record = Entrez.read(handle)

    for datatype in record[0]:
        dictionary[inputID][datatype] = record[0][datatype]

#

print("LOCUS\tNAME\tTaxID\tID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies")

for record in dictionary:

    lineage = {d['Rank']:d['ScientificName'] for d in dictionary[record]["LineageEx"] if d['Rank'] in ['superkingdom','phylum', 'class','order','family','genus','species']}

    kingdom = lineage.get('superkingdom', "unknown")
    phylum = lineage.get('phylum', "unknown")
    taxClass = lineage.get('class', "unknown")
    order = lineage.get('order', "unknown")
    family = lineage.get('family', "unknown")
    genus = lineage.get('genus', "unknown")
    species = lineage.get('species', "unknown")

    print(record,dictionary[record]["ScientificName"],dictionary[record]["TaxId"],dictionary[record]["Id"],kingdom,phylum,taxClass,order,family,genus,species,sep="\t")
