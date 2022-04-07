#!/usr/bin/env python3
"""
IgDiscover computes V/D/J gene usage profiles and discovers novel V genes

- Run IgBLAST in parallel (wrapper inspired by igblastwrp).
- Parse IgBLAST output into a tab-separated table
- Group sequences by barcode
- Plot V gene usage
- Discover new V genes given more than one dataset
"""
import ast
import sys
import logging
import pkgutil
import importlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import matplotlib as mpl
import warnings
import resource

from . import cli as cli_package
from .cli import CommandLineError

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
    if arguments is None:
        arguments = sys.argv[1:]
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    subcommand_name = get_subcommand_name(arguments)
    module = importlib.import_module("." + subcommand_name, cli_package.__name__)

    parser = HelpfulArgumentParser(description=__doc__, prog='igdiscover')
    parser.add_argument('--profile', default=False, action='store_true',
        help='Save profiling information to igdiscover.prof')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--debug', action='store_true', default=False, help='Print debug messages')

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser(
        subcommand_name, help=module.__doc__.split("\n", maxsplit=1)[1], description=module.__doc__
    )
    module.add_arguments(subparser)
    args = parser.parse_args(arguments)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    del args.debug

    do_profiling = args.profile
    del args.profile
    if do_profiling:
        import cProfile as profile
        to_run = lambda: profile.runctx('module.main(args)', globals(), dict(module=module, args=args), filename='igdiscover.prof')
        logger.info('Writing profiling data to igdiscover.prof')
    else:
        to_run = lambda: module.main(args)
    try:
        to_run()
    except CommandLineError as e:
        logger.error(e)
        sys.exit(1)
    if sys.platform == 'linux' and not getattr(module, "do_not_show_cpustats", False):
        rself = resource.getrusage(resource.RUSAGE_SELF)
        rchildren = resource.getrusage(resource.RUSAGE_CHILDREN)
        memory_kb = rself.ru_maxrss + rchildren.ru_maxrss
        cpu_time = rself.ru_utime + rself.ru_stime + rchildren.ru_utime + rchildren.ru_stime
        cpu_time_s = format_duration(cpu_time)
        logger.info('CPU time {}. Maximum memory usage {:.3f} GB'.format(
            cpu_time_s, memory_kb / 1E6))


def get_subcommand_name(arguments) -> str:
    """
    Parse arguments to find out which subcommand was requested.

    This sets up a minimal ArgumentParser with the correct help strings.

    Because help is obtained from a module’s docstring, but importing each module
    makes startup slow, the modules are only parsed with the ast module and
    not fully imported at this stage.

    Return:
        subcommand name
    """
    parser = HelpfulArgumentParser(description=__doc__, prog="igdiscover")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    for module_name, docstring in cli_modules(cli_package):
        help = docstring.split("\n", maxsplit=1)[1].replace("%", "%%")
        subparser = subparsers.add_parser(
            module_name, help=help, description=docstring, add_help=False
        )
        subparser.set_defaults(module_name=module_name)
    args, _ = parser.parse_known_args(arguments)
    module_name = getattr(args, "module_name", None)
    if module_name is None:
        parser.error("Please provide the name of a subcommand to run")
    return module_name


def cli_modules(package):
    """
    Yield (module_name, docstring) tuples for all modules in the given package.
    """
    modules = pkgutil.iter_modules(package.__path__)
    for module in modules:
        spec = importlib.util.find_spec(package.__name__ + "." + module.name)
        with open(spec.origin) as f:
            mod_ast = ast.parse(f.read())
        docstring = ast.get_docstring(mod_ast, clean=False)
        yield module.name, docstring


if __name__ == '__main__':
    main()
