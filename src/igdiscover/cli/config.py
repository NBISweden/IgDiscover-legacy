"""
Change configuration file

Example:
    igdiscover config --set iterations 1

If no --set option is given, print the current configuration.
"""
import os
import logging
import sys

from ruamel.yaml import YAML

from ..config import Config

logger = logging.getLogger(__name__)
do_not_show_cpustats = 1


def add_arguments(parser):
    arg = parser.add_argument
    arg('--set', nargs=2, default=[], metavar=('KEY', 'VALUE'), action='append',
        help='Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys.'
            'Can be given multiple times.')
    arg('--file', default=Config.DEFAULT_PATH,
        help='Configuration file to modify. Default: igdiscover.yaml in current directory.')


def main(args):
    if args.set:
        modify_configuration(args.set, args.file)
    else:
        print_configuration(args.file)


def modify_configuration(
    settings,
    path=Config.DEFAULT_PATH,
):
    yaml = YAML()
    with open(path) as f:
        config = yaml.load(f)
    for k, v in settings:
        if not isinstance(k, str) or not isinstance(v, str):
            raise ValueError("key and value must both be strings")
        v = yaml.load(v)
        # config[k] = v
        item = config
        # allow nested keys
        keys = k.split('.')
        for i in keys[:-1]:
            item = item[i]
        item[keys[-1]] = v
    tmpfile = path + '.tmp'
    with open(tmpfile, 'w') as f:
        yaml.dump(config, f)
    os.rename(tmpfile, path)


def print_configuration(path=Config.DEFAULT_PATH):
    yaml = YAML(typ="safe", pure=True)
    yaml.default_flow_style = False
    with open(path) as f:
        config = yaml.load(f)
    yaml.dump(config, sys.stdout)
