import argparse
import json
import jak_utils

from jakomics import patric

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description='XXX', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with top-level patric folders",
                    required=True)

parser.add_argument('-s', '--sep',
                    help="Column separator",
                    default="\t",
                    required=False)

args = parser.parse_args()

# MAIN ########################################################################

if __name__ == "__main__":

    jak_utils.header()

    files = patric.get_files(args.in_dir)

    print("GENOME", "SUPERCLASS", "CLASS", "SUBCLASS", "SUBSYSTEM",
          "ACTIVE", "PRODUCT", "LOCUS_TAG", sep=args.sep)

    if "|" in args.sep:
        print("---", "---", "---", "---", "---",
              "---", "---", "---", sep=args.sep)

    for file in files:
        subsystem = json.loads(open(file.subsystems_path, "r").read())
        for a in subsystem:
            s = patric.Subsystem(a)
            s.print_hierarchy(sep=args.sep)
