import yaml
import sys
import os
import importlib.metadata
import datetime
from jakomics import colors


def header(r = False):
    code_name = os.path.basename(sys.argv[0])
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    j_version = importlib.metadata.version("jakomics")

    if r == False:
        print(
            f'{colors.bcolors.BLUE}{code_name} - JAK-BIO v{data["version"]} - JAKomics v{j_version}{colors.bcolors.END}', file=sys.stderr)
    else:
        return [code_name,  f'JAK-BIO v{data["version"]}', f'JAKomics v{j_version}', f'START_TIME {timestamp()}']


def get_yaml(field):
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    return data[field]


def timestamp():
    return datetime.datetime.now().strftime("%Y%m%d %H:%M:%S")
