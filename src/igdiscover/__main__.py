#!/usr/bin/env python3
"""
IgDiscover computes V/D/J gene usage profiles and discovers novel V genes

- Run IgBLAST in parallel (wrapper inspired by igblastwrp).
- Parse IgBLAST output into a tab-separated table
- Group sequences by barcode
- Plot V gene usage
- Discover new V genes given more than one dataset
"""
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib as mpl
import warnings
import resource

import igdiscover.cli as cli_package
from igdiscover.cli import CommandLineError

from . import __version__

__author__ = "Marcel Martin"

mpl.use('Agg')
warnings.filterwarnings('ignore', 'axes.color_cycle is deprecated and replaced with axes.prop_cycle')
warnings.filterwarnings('ignore', 'The `IPython.html` package')
warnings.filterwarnings('ignore', 'Widget registration using a string name has been deprecated')
warnings.filterwarnings('ignore', 'Traits should be given as instances, not types')
warnings.filterwarnings('ignore', 'pandas.util.testing is deprecated.')

logger = logging.getLogger(__name__)


class HelpfulArgumentParser(ArgumentParser):
    """An ArgumentParser that prints full help on errors."""

    def __init__(self, *args, **kwargs):
        if 'formatter_class' not in kwargs:
            kwargs['formatter_class'] = RawDescriptionHelpFormatter
        super().__init__(*args, **kwargs)

    def error(self, message):
        self.print_help(sys.stderr)
        args = {'prog': self.prog, 'message': message}
        self.exit(2, '%(prog)s: error: %(message)s\n' % args)


def format_duration(seconds):
    h = int(seconds // 3600)
    seconds -= h * 3600
    m = int(seconds // 60)
    seconds -= m * 60
    return '{:02d}:{:02d}:{:04.1f}'.format(h, m, seconds)


def main(arguments=None):
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    parser = HelpfulArgumentParser(description=__doc__, prog='igdiscover')
    parser.add_argument('--profile', default=False, action='store_true',
        help='Save profiling information to igdiscover.prof')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--debug', action='store_true', default=False, help='Print debug messages')

    show_cpustats = dict()
    subparsers = parser.add_subparsers()

    # Import each module that implements a subcommand and add a subparser for it.
    # Each subcommand is implemented as a module in the cli subpackage.
    # It needs to implement an add_arguments() and a main() function.
    modules = pkgutil.iter_modules(cli_package.__path__)
    for _, module_name, _ in modules:
        module = importlib.import_module("." + module_name, cli_package.__name__)
        subparser = subparsers.add_parser(module_name,
            help=module.__doc__.split('\n')[1], description=module.__doc__)
        subparser.set_defaults(func=module.main)
        module.add_arguments(subparser)
        if hasattr(module, 'do_not_show_cpustats'):
            show_cpustats[module.main] = False

    args = parser.parse_args(arguments)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    do_profiling = args.profile
    del args.profile
    if hasattr(args, 'func'):
        subcommand = args.func
        del args.func
    else:
        parser.error('Please provide the name of a subcommand to run')
    if do_profiling:
        import cProfile as profile
        to_run = lambda: profile.runctx('subcommand(args)', globals(), dict(subcommand=subcommand, args=args), filename='igdiscover.prof')
        logger.info('Writing profiling data to igdiscover.prof')
    else:
        to_run = lambda: subcommand(args)
    try:
        to_run()
    except CommandLineError as e:
        logger.error(e)
        sys.exit(1)
    if sys.platform == 'linux' and show_cpustats.get(subcommand, True):
        rself = resource.getrusage(resource.RUSAGE_SELF)
        rchildren = resource.getrusage(resource.RUSAGE_CHILDREN)
        memory_kb = rself.ru_maxrss + rchildren.ru_maxrss
        cpu_time = rself.ru_utime + rself.ru_stime + rchildren.ru_utime + rchildren.ru_stime
        cpu_time_s = format_duration(cpu_time)
        logger.info('CPU time {}. Maximum memory usage {:.3f} GB'.format(
            cpu_time_s, memory_kb / 1E6))


if __name__ == '__main__':
    main()
