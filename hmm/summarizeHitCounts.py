import sys
import statistics

minScore = 25
geneDict = {}

lines = [line.strip() for line in open(sys.argv[1])]

for line in list(lines):
    if not line.startswith("#"):
        split = line.split()

        if len(split) >= 13:
            model = split[0]
            gene = split[3]
            score = float(split[13])

            if score >= minScore:
                if gene in geneDict:
                    #print(gene,geneDict[gene],sep="\t")
                    if score > geneDict[gene]["score"]:
                        geneDict[gene]["model"] = model
                        geneDict[gene]["score"] = score
                else:
                    geneDict[gene] = {"model" : model, "score" : score}

# now, get the sum of the final hits

modelDict = {}

for gene in geneDict:
    model = geneDict[gene]["model"]
    modelDict[model] = modelDict.get(model, 0) + 1

for model in modelDict:
    print(model,modelDict[model],sep="\t")

print("MEDIAN: "+str(statistics.median(list(modelDict.values()))))
print("MEAN: "+str(statistics.mean(list(modelDict.values()))))
