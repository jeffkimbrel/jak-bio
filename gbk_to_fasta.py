from Bio import SeqIO
import re
import os
import argparse
import sys

from jakomics.gene import GENE

# RAST/patric output by default has an incorrect header which leads to a BioPython warning. This will suppress all warnings.
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def main(gbk, write_contig = None, write_faa = None, write_nt = None, feature_identifier = 'locus_tag', return_gene_dict = False):

    if write_contig != None:
        ct_file = open(write_contig, 'w')

    if write_faa != None:
        aa_file = open(write_faa, 'w')

    if write_nt != None:
        nt_file = open(write_nt, 'w')

    if return_gene_dict:
        gene_dict = {}

    for seq_record in SeqIO.parse(gbk, "genbank"):
        if write_contig != None:
            ct = str(seq_record.seq)
            ct_file.write(">" + seq_record.name + "\n" + ct + "\n")

        for feature in seq_record.features:
            if feature.type in ['CDS', 'ncRNA']:
                id = None
                if feature_identifier in feature.qualifiers:
                    id = feature.qualifiers[feature_identifier][0]
                elif 'db_xref' in feature.qualifiers:
                    for f in feature.qualifiers['db_xref']:
                        if f.startswith('RAST2:fig|'):
                            id = f.replace('RAST2:fig|', '')

                if id is None:
                    print("ERROR: NO FASTA IDENTIFIER FOUND", file=sys.stderr)

                else:
                    if return_gene_dict:
                        gene = GENE(id)
                        gene.replicon = seq_record.name
                        gene.parse_gbk_location(feature.location)

                    if write_faa != None:
                        if 'translation' in feature.qualifiers:
                            aa = feature.qualifiers['translation'][0]
                            aa_file.write(">" + id + "\n" + aa + "\n")

                            if return_gene_dict:
                                gene.aa = aa

                    if write_nt != None:
                        nt = str(feature.location.extract(seq_record).seq)
                        nt_file.write(">" + id + "\n" + nt + "\n")

                        if return_gene_dict:
                            gene.nt = nt

                    if return_gene_dict:
                        gene_dict[id] = gene

    if write_contig != None:
        ct_file.close()

    if write_faa != None:
        aa_file.close()

    if write_nt != None:
        nt_file.close()

    if return_gene_dict:
        return gene_dict

# MAIN ########################################################################

if __name__ == "__main__":

    # OPTIONS #################################################################

    parser = argparse.ArgumentParser(description='Create fasta files from a genbank record')

    parser.add_argument('-g', '--genbank',
                        help="Genbank file to summarize",
                        required=True)

    parser.add_argument('-i', '--identifier',
                        default="locus_tag",
                        help="Which genbank identifier should be the fasta header")

    parser.add_argument('-a', '--amino',
                        default=None,
                        help="Write amino acid fasta file",
                        required=False)

    parser.add_argument('-n', '--nucleo',
                        default=None,
                        help="Write nucleic acid fasta file",
                        required=False)

    parser.add_argument('-c', '--contig',
                        default=None,
                        help="Write contig fasta file",
                        required=False)

    parser.add_argument('--return_gene_dict',
                        default=False,
                        help="Return a list of gene classes",
                        required=False)

    args = parser.parse_args()

    main(args.genbank, args.contig, args.amino, args.nucleo, args.identifier, args.return_gene_dict)
