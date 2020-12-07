import argparse
import pandas as pd
import os
from jakomics import utilities, colors

# FUNCTIONS ###################################################################


def get_counts(file_path, column, header, remove_text):
    if header == 'None':
        column = int(column)
        df = pd.read_csv(file_path, sep="\t", header=None, names=[0, 1])
    else:
        df = pd.read_csv(file_path, sep="\t")

    df = df[column].value_counts()
    df.name = os.path.basename(file_path).replace(remove_text, "")
    return pd.DataFrame(df)


def main(file_list, file_out, column, header, remove_text):
    '''
    Takes a list of table files, and merges on a column
    '''
    shared_df = pd.DataFrame()

    for file in file_list:
        print(f'Merging {column} results from {colors.bcolors.GREEN}{file}{colors.bcolors.END}')

        try:
            shared_df = pd.merge(shared_df, get_counts(
                file, column, header, remove_text), how='outer', left_index=True, right_index=True)
        except IndexError:
            shared_df = shared_df.reindex_axis(
                shared_df.columns.union(get_counts(file, column, header, remove_text).columns), axis=1)

    shared_df = shared_df.fillna(0)
    shared_df.to_csv(file_out, sep="\t")
    print(f'\nSaved merged table to {colors.bcolors.GREEN}{file_out}{colors.bcolors.END}')


# MAIN ########################################################################


if __name__ == "__main__":

    # OPTIONS #####################################################################

    print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

    parser = argparse.ArgumentParser(description="",
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

    parser.add_argument('-t', '--text',
                        help="Common file name text to identify files",
                        nargs='*',
                        required=False,
                        default=['.txt'])

    parser.add_argument('-r', '--remove_text',
                        help="Text to remove from sample names",
                        required=False,
                        default="")

    parser.add_argument('-c', '--column',
                        help="Column name to count",
                        required=True)

    parser.add_argument('-o', '--out',
                        help="Output file name",
                        default='count_and_merge_out.txt')

    parser.add_argument('--header',
                        help='Row number for header. Use None if no header',
                        required=False,
                        default=1)

    args = parser.parse_args()

    file_list = utilities.get_file_list(args.files, ending=args.text)
    file_list = utilities.get_directory_file_list(
        args.in_dir, ending=args.text, file_list=file_list)

    main(file_list, args.out, args.column, args.header, args.remove_text)
