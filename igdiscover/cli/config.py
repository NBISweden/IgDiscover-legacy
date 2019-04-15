"""
Change configuration file

Example:
    igdiscover config --set iterations 1

If no --set option is given, print the current configuration.
"""
import os
import logging
import ruamel.yaml

from igdiscover.config import Config

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--set', nargs=2, default=[], metavar=('KEY', 'VALUE'), action='append',
        help='Set KEY to VALUE. Use KEY.SUBKEY[.SUBSUBKEY...] for nested keys.'
            'Can be given multiple times.')
    arg('--file', default=Config.DEFAULT_PATH,
        help='Configuration file to modify. Default: igdiscover.yaml in current directory.')


def main(args):
    if args.set:
        with open(args.file) as f:
            config = ruamel.yaml.load(f, ruamel.yaml.RoundTripLoader)
        for k, v in args.set:
            v = ruamel.yaml.safe_load(v)
            # config[k] = v
            item = config
            # allow nested keys
            keys = k.split('.')
            for i in keys[:-1]:
                item = item[i]
            item[keys[-1]] = v
        tmpfile = args.file + '.tmp'
        with open(tmpfile, 'w') as f:
            print(ruamel.yaml.dump(config, Dumper=ruamel.yaml.RoundTripDumper), end='', file=f)
        os.rename(tmpfile, args.file)
    else:
        with open(args.file) as f:
            config = ruamel.yaml.safe_load(f)
        print(ruamel.yaml.dump(config), end='')
