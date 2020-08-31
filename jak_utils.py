import yaml
from jakomics import colors

def header():
    with open('meta.yml') as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        print(f'{colors.bcolors.GREEN}JAK-BIO v{data["version"]}{colors.bcolors.END}')
