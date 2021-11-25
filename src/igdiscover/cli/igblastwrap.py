"""
Run IgBLAST and output the AIRR-formatted result table
"""
import sys
import time
from contextlib import ExitStack
from itertools import islice
import errno
import logging

import dnaio

from ..igblast import IgBlastCache, igblast_parallel_chunked
from ..utils import available_cpu_count
from ..config import GlobalConfig


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg(
        "--threads",
        "-t",
        "-j",
        type=int,
        default=1,
        help="Number of threads. Default: 1. Use 0 for no. of available CPUs.",
    )
    arg("--cache", action="store_true", default=None, help="Use the cache")
    arg(
        "--no-cache",
        action="store_false",
        dest="cache",
        default=None,
        help="Do not use the cache",
    )
    arg(
        "--penalty",
        type=int,
        choices=(-1, -2, -3, -4),
        default=None,
        help="BLAST mismatch penalty (default: -1)",
    )
    arg(
        "--species",
        default=None,
        help="Tell IgBLAST which species to use. Note that this setting does "
        "not seem to have any effect since we provide our own database to "
        "IgBLAST. Default: Use IgBLAST’s default",
    )
    arg(
        "--sequence-type",
        default="Ig",
        choices=("Ig", "TCR"),
        help="Sequence type. Default: %(default)s",
    )
    arg("--limit", type=int, metavar="N", help="Limit processing to first N records")
    arg(
        "--rename",
        default=None,
        metavar="PREFIX",
        help="Rename reads to PREFIXseqN (where N is a number starting at 1)",
    )

    arg("database", help="Database directory with V.fasta, D.fasta, J.fasta.")
    arg("fasta", help="File with original reads")


def main(args):
    config = GlobalConfig()
    use_cache = config.use_cache
    if args.cache is not None:
        use_cache = args.cache
    if use_cache:
        from .. import igblast as igblast_module

        igblast_module._igblastcache = IgBlastCache()  # FIXME
        logger.info("IgBLAST cache enabled")
    if args.threads == 0:
        args.threads = available_cpu_count()
    logger.info("Running IgBLAST on input reads")
    start_time = time.time()
    last_status_update = 0
    with ExitStack() as stack:
        sequences = stack.enter_context(dnaio.open(args.fasta))
        sequences = islice(sequences, 0, args.limit)

        n = 0  # number of records processed so far
        for record in igblast_parallel_chunked(
            args.database,
            sequences,
            sequence_type=args.sequence_type,
            species=args.species,
            threads=args.threads,
            penalty=args.penalty,
            use_cache=use_cache,
        ):
            lines = record.splitlines()
            try:
                if n == 0:
                    print(*lines, sep="\n")
                else:
                    print(*lines[1:], sep="\n")
            except IOError as e:
                if e.errno == errno.EPIPE:
                    sys.exit(1)
                raise
            n += len(lines) - 1
            if n % 1000 == 0:
                elapsed = time.time() - start_time
                if elapsed >= last_status_update + 60:
                    logger.info(
                        "Processed {:10,d} sequences at {:.3f} ms/sequence".format(
                            n, elapsed / n * 1e3
                        )
                    )
                    last_status_update = elapsed
    elapsed = time.time() - start_time
    logger.info(
        "Processed {:10,d} sequences at {:.1f} ms/sequence".format(n, elapsed / n * 1e3)
    )
