from jakomics import utilities, blast, hmm, gene, colors
from jakomics.genome import GENOME
import argparse
import os
import re
import sys
import yaml
import pandas as pd
from natsort import natsorted
# from Bio import SearchIO
from collections import Counter
import math
from multiprocessing import Manager, Pool
from tqdm import tqdm

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa genomes",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

# PREP DATABASES ##############################################################

tadb_data = pd.read_excel(jak_utils.get_yaml("TADB2"), sheet_name="merged", engine='openpyxl')
tadb_type = pd.read_excel(jak_utils.get_yaml("TADB2"), sheet_name="type", engine='openpyxl')
tadb_type = tadb_type.set_index('TA_ID')

blast.make_blast_db('prot', jak_utils.get_yaml("TADB2_faa"))
blast.make_blast_db('nucl', jak_utils.get_yaml("TADB2_ffn"))

# CLASSES #####################################################################


class PAIR():

    def __init__(self, gene1, gene2, id):
        self.id = id
        self.gene1 = gene1
        self.gene2 = gene2

    def get_replicon(self):
        return self.gene1.replicon

    def get_range(self):
        smallest = min(self.gene1.start, self.gene1.stop, self.gene2.start, self.gene2.stop)
        largest = max(self.gene1.start, self.gene1.stop, self.gene2.start, self.gene2.stop)

        return(str(smallest) + "-" + str(largest))

    def get_type(self):
        gene1_type = max(set(self.gene1.tadb_type), key=self.gene1.tadb_type.count)
        gene2_type = max(set(self.gene2.tadb_type), key=self.gene2.tadb_type.count)

        if gene1_type == gene2_type:
            return gene1_type
        else:
            return [gene1_type, gene2_type]

    def all_toxin_tabd_ids(self):
        if self.toxin == self.gene1.id:
            return self.gene1.tadb_ids
        elif self.toxin == self.gene2.id:
            return self.gene2.tadb_ids

    def all_antitoxin_tabd_ids(self):
        if self.antitoxin == self.gene1.id:
            return self.gene1.tadb_ids
        elif self.antitoxin == self.gene2.id:
            return self.gene2.tadb_ids

    def get_all_toxin_names(self):
        if self.toxin == self.gene1.id:
            return dict(Counter(self.gene1.tadb_gene_name))
        elif self.toxin == self.gene2.id:
            return dict(Counter(self.gene2.tadb_gene_name))

    def get_all_antitoxin_names(self):
        if self.antitoxin == self.gene1.id:
            return dict(Counter(self.gene1.tadb_gene_name))
        elif self.antitoxin == self.gene2.id:
            return dict(Counter(self.gene2.tadb_gene_name))

    def view(self, genome_name):
        print(genome_name,
              self.id,
              self.get_replicon(),
              self.get_range(),
              self.get_type(),
              self.toxin,
              self.toxin_name,
              self.antitoxin,
              self.antitoxin_name,
              self.shared_tadb_families,
              self.shared_tadb_ids,
              self.score,
              scores[self.score],
              self.all_toxin_tabd_ids(),
              self.all_antitoxin_tabd_ids(),
              self.get_all_toxin_names(),
              self.get_all_antitoxin_names(),
              sep="\t"
              )

# FUNCTIONS ###################################################################


def gene_distance(gene1, gene2):

    if gene1.id >= gene2.id:  # only compare once
        return math.inf
    elif gene1.replicon != gene2.replicon:
        return math.inf
    else:
        shortest_distance = min(
            abs(gene1.start - gene2.start),
            abs(gene1.stop - gene2.stop),
            abs(gene1.stop - gene2.start),
            abs(gene1.start - gene2.stop)
        )
        return shortest_distance


def get_tadb_family(tadb):
    tadb = tadb.replace("TADB|", "")
    role, ta_id, _ = re.split('(\d+)', tadb)
    df = tadb_data.loc[tadb_data['TA_ID'] == int(ta_id)]
    family = df['TA_FAMILY'].tolist()
    prediction = df[role].tolist()

    return family, prediction, role, ta_id


def score_tadb_pairs(pairs):
    # print(f'Scoring Potential Pairs')
    processed_pairs = []

    for pair in pairs:
        pair.gene1_type = max(set(pair.gene1.tadb_type), key=pair.gene1.tadb_type.count)
        pair.gene1_role = max(set(pair.gene1.tadb_roles), key=pair.gene1.tadb_roles.count)
        pair.gene1_family = max(set(pair.gene1.tadb_families), key=pair.gene1.tadb_families.count)
        pair.gene1_name = max(set(pair.gene1.tadb_gene_name), key=pair.gene1.tadb_gene_name.count)
        pair.shared_tadb_ids = {}
        pair.shared_tadb_families = {}

        if pair.gene2 != None:
            pair.gene2_type = max(set(pair.gene2.tadb_type), key=pair.gene2.tadb_type.count)
            pair.gene2_role = max(set(pair.gene2.tadb_roles), key=pair.gene2.tadb_roles.count)
            pair.gene2_family = max(set(pair.gene2.tadb_families),
                                    key=pair.gene2.tadb_families.count)
            pair.gene2_name = max(set(pair.gene2.tadb_gene_name),
                                  key=pair.gene2.tadb_gene_name.count)
            pair.shared_tadb_ids = list(
                set(pair.gene1.tadb_ids).intersection(set(pair.gene2.tadb_ids)))
            pair.shared_tadb_families = list(
                set(pair.gene1.tadb_families).intersection(set(pair.gene2.tadb_families)))

            if pair.gene1_role == "T":
                pair.toxin = pair.gene1.id
                pair.toxin_name = pair.gene1_name
                pair.antitoxin = pair.gene2.id
                pair.antitoxin_name = pair.gene2_name
            else:
                pair.toxin = pair.gene2.id
                pair.toxin_name = pair.gene2_name
                pair.antitoxin = pair.gene1.id
                pair.antitoxin_name = pair.gene1_name

            if pair.gene1_role == pair.gene2_role:
                if len(set(pair.gene1.tadb_roles)) > 1 or len(set(pair.gene2.tadb_roles)) > 1:
                    pair.gene1_role = dict(Counter(pair.gene1.tadb_roles))
                    pair.gene2_role = dict(Counter(pair.gene2.tadb_roles))
                    pair.score = 5
                else:
                    pair.score = 1

            else:
                # start more strict
                if len(pair.shared_tadb_ids) > 0:
                    pair.score = 20
                elif pair.gene1_family == pair.gene2_family:
                    pair.score = 15
                elif len(pair.shared_tadb_families) > 0:
                    pair.score = 10
                else:
                    pair.score = 0

        processed_pairs.append(pair)

    return(pairs)


scores = {1:  "Not candidates - have same role",
          0:  "Same TYPE and opposite ROLEs",
          5:  "Some ambiguity in ROLEs",
          10: "Have at least one TADB FAMILY prediction in common",
          15: "Best TADB FAMILY predictions are the same",
          20: "TA predictions have TADB IDs in common"}


def parse_tadb_results(gene):
    gene.tadb_families = []
    gene.tadb_gene_name = []
    gene.tadb_roles = []
    gene.tadb_ids = []
    gene.tadb_type = []

    if hasattr(gene, 'tadb_blast'):
        for blast_result in gene.tadb_blast:
            if blast_result.filter(e=1e-15, p=25):
                family, prediction, role, ta_id = get_tadb_family(blast_result.subject)
                gene.tadb_families += family
                gene.tadb_gene_name += prediction
                gene.tadb_roles.append(role)
                gene.tadb_ids.append(ta_id)
                gene.tadb_type.append(tadb_type.loc[int(ta_id), 'TA_TYPE'])

    return gene


def get_potential_tadb_pairs(genome, overlapping_ids=1, max_distance=500):
    # print(f'Finding Potential Pairs')

    pairs = []

    genome.potential_TA_list = list(set(genome.potential_TA_list))
    for gene1 in natsorted(genome.potential_TA_list):
        for gene2 in natsorted(genome.potential_TA_list):
            distance = gene_distance(genome.genes[gene1], genome.genes[gene2])
            if distance <= max_distance:
                pairs.append(PAIR(genome.genes[gene1], genome.genes[gene2], len(pairs) + 1))

     # add genes without pairs here

    return pairs


def blast_tadb(genome, aa=True, nt=False):
    # print(f'Starting BLASTs')

    if aa:
        tadb_aa_results = blast.run_blast(type="prot",
                                          q=genome.faa_path,
                                          db=jak_utils.get_yaml("TADB2_faa"),
                                          e=1e-7,
                                          make=False)

        for query in natsorted(tadb_aa_results.keys()):
            genome.genes[query].tadb_blast = tadb_aa_results[query]
            # for hit in tadb_aa_results[query]:
            #     hit.print_rough_result()

    if nt:
        tadb_nt_results = blast.run_blast(type="nucl",
                                          q=genome.nt_path,
                                          db=jak_utils.get_yaml("TADB2_ffn"),
                                          e=1e-2,
                                          make=False)

        for query in natsorted(tadb_nt_results.keys()):
            genome.genes[query].tadb_blast = tadb_nt_results[query]
            # for hit in tadb_nt_results[query]:
            #     hit.print_rough_result()
    #
    #
    # tadb_contig_results = blast.run_blast(type="nucl",
    #                           q=genome.contig_path,
    #                           db="/Users/kimbrel1/Science/repos/jakomics/db/TADB_2.0/TADB2.ffn",
    #                           e=1e-7,
    #                           make = True)
    #
    # for query in natsorted(tadb_contig_results.keys()):
    #     genome.genes[query].tadb.append(tadb_contig_results[query])


def view_tadb_results(name, pairs, min_score):
    print('GENOME', 'GENOME_PAIR', 'REPLICON', 'RANGE', 'TYPE', 'TOXIN_LOCUS', 'TOXIN_NAME', 'ANTITOXIN_LOCUS', 'ANTITOXIN_NAME', 'SHARED_FAMILIES',
          'SHARED_TADB_IDS', 'SCORE', 'SCORE_NOTE', 'TOXIN_TADB_IDS', 'ANTITOXIN_TADB_IDS', 'TOXIN_NAMES', 'ANTITOXIN_NAMES', sep="\t")

    for pair in pairs:
        if pair.score >= min_score:
            pair.view(name)


def make_empty_df():
    return(pd.DataFrame(columns=['GENOME', 'GENOME_PAIR', 'REPLICON', 'RANGE', 'TYPE', 'TOXIN_LOCUS', 'TOXIN_NAME', 'ANTITOXIN_LOCUS', 'ANTITOXIN_NAME', 'SHARED_FAMILIES',
                                 'SHARED_TADB_IDS', 'SCORE', 'SCORE_NOTE', 'TOXIN_TADB_IDS', 'ANTITOXIN_TADB_IDS', 'TOXIN_NAMES', 'ANTITOXIN_NAMES']))


def tadb_results_to_df(name, pairs, min_score):
    details = make_empty_df()

    for pair in pairs:
        if pair.score >= min_score:
            # pair.view(name)
            details = details.append(
                pd.Series(data={
                    'GENOME': name,
                    'GENOME_PAIR': pair.id,
                    'REPLICON': pair.get_replicon(),
                    'RANGE': pair.get_range(),
                    'TYPE': pair.get_type(),
                    'TOXIN_LOCUS': pair.toxin,
                    'TOXIN_NAME': pair.toxin_name,
                    'ANTITOXIN_LOCUS': pair.antitoxin,
                    'ANTITOXIN_NAME': pair.antitoxin_name,
                    'SHARED_FAMILIES': pair.shared_tadb_families,
                    'SHARED_TADB_IDS': pair.shared_tadb_ids,
                    'SCORE': pair.score,
                    'SCORE_NOTE': scores[pair.score],
                    'TOXIN_TADB_IDS': pair.all_toxin_tabd_ids(),
                    'ANTITOXIN_TADB_IDS': pair.all_antitoxin_tabd_ids(),
                    'TOXIN_NAMES': pair.get_all_toxin_names(),
                    'ANTITOXIN_NAMES': pair.get_all_antitoxin_names()
                }
                ),
                ignore_index=True)

    return details


def find_TAs(genome):

    global shared_results
    results = []

    # write genes to genomes and gene class dictionary
    gbk = GENOME(genome)

    # write genes to genomes and gene class dictionary
    genome.faa_path = os.path.join(args.out_dir, genome.name + ".faa")
    genome.nt_path = os.path.join(args.out_dir, genome.name + ".ffn")
    genome.contig_path = os.path.join(args.out_dir, genome.name + ".fa")

    genome.genes = gbk.genbank_to_fasta(write_faa=genome.faa_path,
                                        write_nt=genome.nt_path,
                                        write_contig=genome.contig_path,
                                        return_gene_dict=True)

    genome.potential_TA_list = []

    blast_tadb(genome, aa=True, nt=True)

    for gene in natsorted(genome.genes):
        genome.genes[gene] = parse_tadb_results(genome.genes[gene])

        if len(genome.genes[gene].tadb_roles) > 0:
            genome.potential_TA_list.append(genome.genes[gene].id)

    pairs = get_potential_tadb_pairs(genome)
    shared_results[genome.name] = score_tadb_pairs(pairs)

# MAIN LOOP ###################################################################


if __name__ == "__main__":
    jak_utils.header()

    manager = Manager()
    shared_results = manager.dict()

    if not os.path.exists(args.out_dir):
        print("\nCreating directory " + args.out_dir)
        os.makedirs(args.out_dir)

    genome_list = utilities.get_files(args.files, args.in_dir, ["gbk", "gbff", "gb"])

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(find_TAs, genome_list), total=len(genome_list), desc="Finished", unit=" genomes"):
        pass
    pool.close()

    results_df = make_empty_df()
    for genome in shared_results:
        details = tadb_results_to_df(genome, shared_results[genome], min_score=5)
        results_df = pd.concat([results_df, details])

    print(results_df)

    # write to file with comments
    f = open(os.path.join(args.out_dir, "results.txt"), 'a')
    for c in jak_utils.header(r=True):
        print(f'# {c}', file=f)
    for arg in vars(args):
        print(f'# ARG {arg} = {getattr(args, arg)}', file=f)

    results_df.to_csv(f, sep="\t", index=False)
