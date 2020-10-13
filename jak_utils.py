import yaml
import sys
import os
from jakomics import colors


def header():
    code_name = os.path.basename(sys.argv[0])
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    print(
        f'{colors.bcolors.BLUE}JAK-BIO v{data["version"]}: {code_name}{colors.bcolors.END}', file=sys.stderr)


def get_yaml(field):
    yaml_path = os.path.join(os.path.dirname(sys.argv[0]), 'jak_meta.yml')
    with open(yaml_path) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)

    return data[field]
