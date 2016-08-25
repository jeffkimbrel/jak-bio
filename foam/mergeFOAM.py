import glob

S = {} # key = sample, value = KEGG:counts
K = {} # key = kegg, value = kegg

for f in glob.glob('*.KO'):
    
    S[f] = {}
    lines = [line.strip() for line in open(f)]
        
    for line in lines:
        kegg,count = line.split()
        K[kegg] = kegg
        S[f][kegg] = count

sortedK = sorted(list(K))
sortedS = sorted(list(S))

# print header

for sample in sortedS:
    split = sample.split("_")
    print("\t"+split[0],end="")

print()


for kegg in sortedK:
    print(kegg,end="")
    for sample in sortedS:
        #print(S[sample][kegg])
        print("\t"+str(S[sample].get(kegg, 0)),end="")
    print()
        
