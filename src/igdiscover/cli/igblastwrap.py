"""
Run IgBLAST and output the AIRR-formatted result table
"""
import sys
import os
import time
import multiprocessing
from contextlib import ExitStack
from itertools import islice
import errno
import logging
import tempfile

import dnaio

from ..igblast import chunked, makeblastdb, run_igblast, IgBlastCache
from ..utils import SerialPool, available_cpu_count
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


class Runner:
    """
    This is the target of a multiprocessing pool. The target needs to
    be pickleable, and because nested functions cannot be pickled,
    we need this separate class.

    It runs IgBLAST and parses the output for a list of sequences.
    """

    def __init__(
        self, blastdb_dir, species, sequence_type, penalty, database, use_cache
    ):
        self.blastdb_dir = blastdb_dir
        self.species = species
        self.sequence_type = sequence_type
        self.penalty = penalty
        self.database = database
        self.use_cache = use_cache

    def __call__(self, sequences):
        """
        Return raw IgBLAST output
        """
        return run_igblast(
            sequences,
            self.blastdb_dir,
            self.species,
            self.sequence_type,
            self.penalty,
            self.use_cache,
            airr=True,
        )


def igblast(
    database,
    sequences,
    sequence_type,
    species=None,
    threads=1,
    penalty=None,
    use_cache=False,
):
    """
    Run IgBLAST on chunks of the input and yield the AIRR-formatted results.

    database -- Path to database directory with V./D./J.fasta files
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    threads -- number of threads.
    """
    with ExitStack() as stack:
        # Create the three BLAST databases in a temporary directory
        blastdb_dir = stack.enter_context(tempfile.TemporaryDirectory())
        for gene in ["V", "D", "J"]:
            # Without adding the "%" prefix, IgBLAST reports record names that look like GenBank
            # ids as "gb|original_name|".
            makeblastdb(
                os.path.join(database, gene + ".fasta"),
                os.path.join(blastdb_dir, gene),
                prefix="%",
            )

        chunks = chunked(sequences, chunksize=1000)
        runner = Runner(
            blastdb_dir,
            species=species,
            sequence_type=sequence_type,
            penalty=penalty,
            database=database,
            use_cache=use_cache,
        )
        pool = stack.enter_context(
            multiprocessing.Pool(threads) if threads > 1 else SerialPool()
        )
        for igblast_output in pool.imap(runner, chunks, chunksize=1):
            yield igblast_output


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
        for record in igblast(
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
