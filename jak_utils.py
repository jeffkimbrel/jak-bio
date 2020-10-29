import yaml
import sys
import os
import pkg_resources

from jakomics import colors


def header(r=False):
    code_name = os.path.basename(sys.argv[0])
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    j_version = pkg_resources.get_distribution("jakomics").version

    if r == False:
        print(
            f'{colors.bcolors.BLUE}{code_name} - JAK-BIO v{data["version"]} - JAKomics v{j_version}{colors.bcolors.END}', file=sys.stderr)
    else:
        return [code_name,  f'JAK-BIO v{data["version"]}', f'JAKomics v{j_version}']


def get_yaml(field):
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    return data[field]
