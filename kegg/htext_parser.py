import sys
import re

lines = [line.strip() for line in open(sys.argv[1])]

hier = {}

for line in list(lines):
    if re.match("^[A-E]", line):

        if line.startswith("A"):
            A = line[4:-4]
            if len(A) > 0:
                hier["A"] = A

        elif line.startswith("B"):
            B = line[1:]

            B = B.replace("<b>", "")
            B = B.replace("</b>", "")
            B = B.lstrip()
            B = B.rstrip()

            if len(B) > 0:
                hier["B"] = B

        elif line.startswith("C"):
            C = line[1:]

            C = C.replace("<b>", "")
            C = C.replace("</b>", "")
            C = C.lstrip()
            C = C.rstrip()

            if len(C) > 0:
                hier["C"] = C

        elif line.startswith("D"):
            D = line[1:]

            D = D.replace("<b>", "")
            D = D.replace("</b>", "")
            D = D.lstrip()
            D = D.rstrip()

            if len(D) > 0:
                module = D[:6]
                pathway = D[8:]
                hier["D"] = {"module" : module, "pathway" : pathway}


        elif line.startswith("E"):
            E = line[1:]

            E = E.replace("<b>", "")
            E = E.replace("</b>", "")
            E = E.lstrip()
            E = E.rstrip()

            if len(E) > 0:
                ko = E[:6]
                product = E[8:]
                keggID = ko+"_"+hier["D"]["module"]
                hier["E"] = {"ko" : ko, "product" : product}

                print(keggID, ko, product, hier["A"], hier["B"], hier["C"], hier["D"]["module"], hier["D"]["pathway"], sep = "\t")
