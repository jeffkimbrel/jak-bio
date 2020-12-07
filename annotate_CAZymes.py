import sys
import os
import argparse
import uuid
from multiprocessing import Manager, Pool
from tqdm import tqdm

from jakomics import hmm, utilities, colors
import count_and_merge_tables

import jak_utils

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='''

            Find CAZYmes in an amino acid file

            QC_CODE refers to different scoring rules set forth by DBCAN6:
                1A = E-value <= 1e-18 AND HMM coverage >= 35%
                1B = E-value <= 1e-15 AND HMM coverage >= 35%
                2 = E-value <= 1e-5  AND Alignment Length >= 80
                3 = E-value <= 1e-3  AND HMM coverage >= 30%
                0 = Low-quality hit that did not pass any metrics''',

                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--qc',
                    help="QC codes to include",
                    required=False,
                    nargs='*',
                    default=['1A', '1B'])

parser.add_argument('--hmm',
                    help="Write merged HMM results",
                    required=False,
                    default=None)

parser.add_argument('--substrate',
                    help="Write merged substrates results",
                    required=False,
                    default=None)

args = parser.parse_args()

manager = Manager()
counter = manager.dict()

# CLASSES #####################################################################

genes = {}


class File:

    def __init__(self, file_path):
        self.file_path = file_path
        self.run_id = utilities.get_unique_ID()
        self.file_name = os.path.basename(file_path)
        self.results_file = self.file_name + '.dbcan8.txt'
        self.results_class_file = self.file_name + '.class_summary.txt'
        # self.results_berlemont_file = self.file_name + '.berlemont_summary.txt'
        self.temp_log = self.run_id + '.log'
        self.temp_output = self.run_id + '.temp.txt'

    def remove_temp(self):
        os.system('rm ' + self.temp_log)
        os.system('rm ' + self.temp_output)

    def process_dbcan_results(self):
        '''
        read the raw HMMER file and produce HMM classes
        '''

        rawResult = [line.strip() for line in open(self.temp_output)]

        self.cazymes = {}

        for line in rawResult:
            if not line.startswith("#"):

                cazyme = hmm.CAZYME(line)
                cazyme.pass_cazy()
                if cazyme.gene in self.cazymes:
                    self.cazymes[cazyme.gene].append(cazyme)
                else:
                    self.cazymes[cazyme.gene] = [cazyme]

    def write_results(self):
        hmm_output = open(self.results_file, 'w')
        hmm_output.write(
            "LOCUS\tHMM\tEVAL\tSCORE\tC-EVAL\tI-EVAL\tSEQ_COORDS\tHMM_COORDS\tALIGN_LENGTH\tHMM_COVERAGE\tQC_CODE\tCLASS\tSUBSTRATE\n")

        for gene in sorted(self.cazymes.keys()):
            for hit in self.cazymes[gene]:
                if hit.pass_qc in args.qc:
                    hit.assign_cazy_class()
                    hit.assign_substrate()
                    hmm_output.write(hit.write())

        hmm_output.close()


# FUNCTIONS ###################################################################

def processHMMER(file):
    '''
    read the raw HMMER file and produce HMM classes
    '''

    rawResult = [line.strip() for line in open(file.temp_output)]

    formattedResult = {}

    for line in rawResult:
        if not line.startswith("#"):

            formattedResult[uuid.uuid4().hex] = hmm.CAZYME(line)

    return(formattedResult)


def main(file_path):

    global counter

    file = File(file_path)

    hmm.run_hmmsearch(file.file_path, file.temp_log, file.temp_output,
                      jak_utils.get_yaml("cazyme_db"))

    file.process_dbcan_results()
    file.write_results()

    # cleanup
    file.remove_temp()
    counter[file.results_file] = file.file_name


## MAIN LOOP ###################################################################


if __name__ == "__main__":
    jak_utils.header()

    file_list = utilities.get_file_list(args.files, [''])
    file_list = utilities.get_directory_file_list(args.in_dir, [''], file_list)

    if len(file_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid .faa files were found!{colors.bcolors.END}")

    pool = Pool(processes=8)
    for _ in tqdm(pool.imap_unordered(main, file_list), total=len(file_list), desc="Finished", unit=" files"):
        pass
    pool.close()

    # write merged results
    if args.hmm is not None:
        count_and_merge_tables.main(counter, args.hmm, 'HMM', 1, '.faa.dbcan8.txt')

    if args.substrate is not None:
        count_and_merge_tables.main(counter, args.substrate, 'HMM', 1, '.faa.dbcan8.txt')
