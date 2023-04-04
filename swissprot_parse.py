from jakomics import annotations
from Bio import SwissProt

# https://biopython.org/docs/1.78/api/Bio.SwissProt.html

import gzip

entry_name_count = 0
accessions_count = 0

with gzip.open("/Users/kimbrel1/OD/Biofuels_SFA/Computational/ProbAnnotation/swiss_prot/uniprot_sprot_2022_01.dat.gz") as handle:
    records = SwissProt.parse(handle)
    for record in records:
        eco = []
        for t in record.comments:
            if t.startswith("CATALYTIC ACTIVITY"):
                if "ECO:0000269" in t:
                    eco.append(t)

        if len(eco) > 0:
            
            for a in record.accessions:
                ecs = []
                ecos = []
                for r in eco:
                    ecs += annotations.find_ec(r)
                    ecos += annotations.parse_evidence_fields(r)

                ecs = list(set(ecs))
                ecos = list(set(ecos))

                ecos = [eco for eco in ecos if eco.startswith("ECO:0000269")]

                if len(ecs) > 0:
                    for ec in ecs:
                        ec = ec.replace("=", ":")
                        # print(f"{a}\t{ec}\t{ecos}")
                        print(f"{a}\t{ec}")
