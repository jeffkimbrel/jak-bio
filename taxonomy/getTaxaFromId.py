from Bio import Entrez
Entrez.email = 'kimbrel1@llnl.gov'

dictionary = {}
inputID = "WP_054761838.1"

dictionary[inputID] = {}


def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = str(taxid), db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

#
handle = Entrez.esearch(db="protein", term=inputID, retmode="xml")
record = Entrez.read(handle)
accession = record["IdList"][0]

#
handle = Entrez.esummary(db="protein", id=accession, retmode="xml")
record = Entrez.read(handle)

tax = record[0]['TaxId']

for item in record[0]:
    dictionary[inputID][item] = record[0][item]

#
handle = Entrez.efetch(db="taxonomy", id=str(tax), retmode="xml")
record = Entrez.read(handle)

for datatype in record[0]:
    dictionary[inputID][datatype] = record[0][datatype]
    
#    
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
    
    
    
    
    #for item in dictionary[record]:
        #print(item,dictionary[record][item],sep=": ")