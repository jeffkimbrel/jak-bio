from pyx import *
import sys
unit.set(defaultunit="mm")
text.set(text.LatexRunner)
text.preamble(r"\usepackage{times}")
c = canvas.canvas()

##### Global Options #####

mmWidth = 1500
mmPadLeft = 30
mmPadRight = 10
mmSpaceHeight = 50
mmFeatureHeight = 20
mmCurrentY = 20

##### Classes #####

class chrome:
    chromeList = []
    def __init__(self, name, start, stop, mmY):
        chrome.chromeList.append(self)
        self.name = name
        self.start = min(start,stop)
        self.stop = max(start,stop)
        self.difference = self.getDiff()
        self.mmY = mmY

    def getDiff(self):
        diff = self.stop - self.start
        return(diff)

class feature:
    featList = []
    def __init__(self, chromosome,start,stop,strand,name):
        feature.featList.append(self)
        self.chromosome = chromosome
        self.start = start
        self.stop = stop
        self.strand = strand
        self.name = name

##### Functions #####

def addChrome(name,start,stop):
    alreadyExists = False
    for instance in chrome.chromeList:
        if instance.name == name:
            instance.start = min(instance.start,start,stop)
            instance.stop = max(instance.stop,start,stop)
            alreadyExists = True
        instance.difference = instance.getDiff()
    if alreadyExists == False:
        chrome(name,start,stop,mmCurrentY)
        return(True) # only to increment the mmCurrentY

def getScaleFactor():
    largestChromeLength = 0
    for instance in chrome.chromeList:
        if instance.difference > largestChromeLength:
            largestChromeLength = instance.difference

    scaleFactor = (mmWidth - mmPadLeft - mmPadRight) / largestChromeLength
    return(scaleFactor)

def drawBackbones():
    largestY = 0
    for instance in chrome.chromeList:
        largestY = max(largestY,instance.mmY)
        start = mmPadLeft-5
        stop = (instance.difference * getScaleFactor()) + start + 10

        c.stroke(path.line(start,instance.mmY, stop, instance.mmY),[style.linewidth(1)])
        c.text(3,instance.mmY-2,instance.name,[text.size.huge])

    boundingBox = path.rect(0, 0, mmWidth, largestY+15)
    c.stroke(boundingBox)

def getChromeStart(chromeName):
    for instance in chrome.chromeList:
        if instance.name == chromeName:
            return(instance.start,instance.mmY)

def drawFeatures():
    for gene in feature.featList:
        chromeStart,mmY = getChromeStart(gene.chromosome)
        start = int((gene.start - chromeStart) * getScaleFactor() + mmPadLeft)
        stop = int((gene.stop - chromeStart) * getScaleFactor() + mmPadLeft) - start
        if gene.strand == '+':
            mmY += (mmFeatureHeight/2)
        elif gene.strand == '-':
            mmY -= (mmFeatureHeight/2)

        c.fill(path.rect(start, mmY-(mmFeatureHeight/2), stop, mmFeatureHeight), [color.cmyk(0.8,0,0,0.4)])

        # I don't know how to print underscores in TeX
        name = gene.name.replace('_',' ')

        c.text(start+3,mmY-2,name,[text.size.normal,color.cmyk.White])

##### Process the file #####

data = open(sys.argv[1], 'rt')
while True:
    line = data.readline()
    if len(line) > 0:
        line = line.rstrip()
        split = line.split('\t')

        incrementY = addChrome(split[0],int(split[1]),int(split[2]))
        if incrementY:
            mmCurrentY += mmSpaceHeight
        feature(split[0],int(split[1]),int(split[2]),split[3],split[4])

    if not line:
        break

##### Actually Draw the Objects #####
drawBackbones()
drawFeatures()

c.writePDFfile("test2.pdf")
