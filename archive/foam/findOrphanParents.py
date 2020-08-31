# this script will verify that the 199 orphan KEGGs are always found with another parent KEGG that is found in the ontology

import sys
import glob

ontology = {}
ontFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/FOAM/FOAM-onto_rel1.tsv")]
for line in ontFile:
    L1,L2,L3,L4,KO = line.split("\t")
    ontology[KO] = KO
    

orphanFile = [line.strip() for line in open("/Users/jak/Dropbox/scripts/FOAM/orphanKOs.txt")]

orphans = {}

for line in orphanFile:
    orphans[line] = {}
    
for f in glob.glob('*.BH'):
    lines = [line.strip() for line in open(f)]
        
    for line in lines:
        split = line.split("\t")
        
        for orphan in orphans:
            if orphan in split[3]:
                koSplit = split[3].split(",")
                
                if len(koSplit) == 1:
                    print(orphan+" is the only KEGG for "+split[4]+"!!!")
                
                for ko in koSplit:
                    
                    junk,kegg = ko.split(":")
                    
                    if kegg in orphans[orphan]:
                        orphans[orphan][kegg] = orphans[orphan][kegg] + 1
                    else:
                        orphans[orphan][kegg] = 1
                        
noParent = {}                   
                        
for orphan in orphans:
    noParent[orphan] = True
    for parent in orphans[orphan]:
        if parent != orphan:
            if parent in ontology:
                print(orphan, str(orphans[orphan][orphan]),parent,str(orphans[orphan][parent]),"YES",sep="\t")
                noParent[orphan] = False
            else:
                print(orphan, str(orphans[orphan][orphan]),parent,str(orphans[orphan][parent]),"NO",sep="\t")
 
 
for orphan in noParent:
    if noParent[orphan] == True:
        print(orphan+" has no parents in the ontology!!!")