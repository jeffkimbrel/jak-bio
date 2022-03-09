import os
import argparse
import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXXXX')

parser.add_argument('-g', '--gtdb',
                    help="gtdb classify output",
                    required=True)

args = parser.parse_args()

# FUNCTIONS ###################################################################

def get_ncbi_dmp():
    print("getting NCBI")
    ncbi = {}
    ncbi_raw = [line.strip() for line in open(os.path.join(jak_utils.get_yaml("ncbi_taxdmp_path"), 'fullnamelineage.dmp'))]
    for line in ncbi_raw:
        split = line.split("|")
        id = int(split[0].strip())

        ncbi[id] = split[2].strip()

    return ncbi

def get_mapping_file():
    print("getting GTDB metadata")
    mapping = {}
    arc = [line.strip() for line in open(os.path.join(jak_utils.get_yaml("gtdb_metadata_path"), 'ar122_metadata_r202.tsv'))]
    bac = [line.strip() for line in open(os.path.join(jak_utils.get_yaml("gtdb_metadata_path"), 'bac120_metadata_r202.tsv'))]
    r202 = arc + bac

    for line in r202:

        split = line.split("\t")
        lineage = split[16]
        taxa = lineage.split(";")[-1]
        id = split[73]  # use 77 for ncbi_taxid, 73 for ncbi species id

        if lineage not in mapping:
            mapping[lineage] = id
        else:
            if id != mapping[lineage] and id != 'none':
                if mapping[lineage] == 'none':
                    mapping[lineage] = id
                else:
                    if int(id) < int(mapping[lineage]):
                        mapping[lineage] = id

    return(mapping)


def get_name(lineage):
    l_split = lineage.split(";")
    for l in l_split.copy():
        if l.endswith('__'):
            # print(l)
            l_split.remove(l)
    return(l_split[-1].split('__')[1])


def find_exact_match(lineage):
    if lineage in gtdb_id:
        return int(gtdb_id[lineage])


def find_partial_match(lineage, gtdb_id):
    # print(lineage)
    l_split = lineage.split(";")
    for l in l_split.copy():
        if l.endswith('__'):
            # print(l)
            l_split.remove(l)
    s = ";".join(l_split)
    # print(s)

    lowest_id = 9999999999

    for id in gtdb_id:
        if s in id:

            if gtdb_id[id] != 'none' and int(gtdb_id[id]) < lowest_id:
                lowest_id = int(gtdb_id[id])

    if lowest_id == 9999999999:
        found = False
        while found == False:
            l_split.pop()
            s = ";".join(l_split)
            for id in gtdb_id:
                if s in id:
                    if gtdb_id[id] != 'none' and int(gtdb_id[id]) < lowest_id:
                        lowest_id = int(gtdb_id[id])
                        found = True

    return int(lowest_id)

# Main ########################################################################


if __name__ == "__main__":
    jak_utils.header()
    gtdb_id = get_mapping_file()
    ncbi = get_ncbi_dmp()

    #
    gtdb = [line.strip() for line in open(args.gtdb)]

    print("GENOME", "NCBI_ID", "NAME", "GTDB_LINEAGE", "NCBI_LINEAGE", sep="\t")

    for line in gtdb:
        split = line.split("\t")
        genome = split[0]
        lineage = split[1]
        # print(line)
        #print(genome + " *" + lineage + "*")

        if genome != 'user_genome':
            tax_id = ''
            if find_exact_match(lineage) != None:
                tax_id = find_exact_match(lineage)
            elif find_partial_match(lineage, gtdb_id) != 9999999999:
                tax_id = find_partial_match(lineage, gtdb_id)
            else:
                print("?????")
            name = get_name(lineage)

            if tax_id in ncbi:
                print(genome, tax_id, name, lineage, ncbi[tax_id], sep="\t")
