"""
Run IgBLAST and output a result table

This is a wrapper for the "igblastn" tool which has a simpler command-line
syntax and can also run IgBLAST in parallel.

The results are parsed, postprocessed and printed as a tab-separated table
to standard output.

Postprocessing includes:

- The CDR3 is detected by using a regular expression
- The leader is detected within the sequence before the found V gene (by
  searching for the start codon).
- If the V sequence hit starts not at base 1 in the reference, it is extended
  to the left.
"""
import sys
import os
import shutil
import time
import multiprocessing
import subprocess
from contextlib import ExitStack
from io import StringIO
from itertools import islice
import hashlib
import errno
import pkg_resources
import logging
import tempfile
import json
import gzip

import dnaio
from xopen import xopen

from ..utils import SerialPool, available_cpu_count, nt_to_aa
from ..parse import TableWriter, IgBlastParser
from ..species import cdr3_start, cdr3_end
from ..config import GlobalConfig


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--threads', '-t', '-j', type=int, default=1,
        help='Number of threads. Default: 1. Use 0 for no. of available CPUs.')
    arg('--cache', action='store_true', default=None, help='Use the cache')
    arg('--no-cache', action='store_false', dest='cache', default=None, help='Do not use the cache')
    arg('--penalty', type=int, choices=(-1, -2, -3, -4), default=None,
        help='BLAST mismatch penalty (default: -1)')
    arg('--species', default=None,
        help='Tell IgBLAST which species to use. Note that this setting does '
            'not seem to have any effect since we provide our own database to '
            'IgBLAST. Default: Use IgBLAST’s default')
    arg('--sequence-type', default='Ig', choices=('Ig', 'TCR'),
        help='Sequence type. Default: %(default)s')
    arg('--raw', metavar='FILE', help='Write raw IgBLAST output to FILE '
            '(add .gz to compress)')
    arg('--limit', type=int, metavar='N',
        help='Limit processing to first N records')
    arg('--rename', default=None, metavar='PREFIX',
        help='Rename reads to PREFIXseqN (where N is a number starting at 1)')
    arg('--stats', metavar='FILE',
        help='Write statistics in JSON format to FILE')

    arg('database', help='Database directory with V.fasta, D.fasta, J.fasta.')
    arg('fasta', help='File with original reads')


class IgBlastCache:
    """Cache IgBLAST results in ~/.cache"""

    binary = 'igblastn'

    def __init__(self):
        version_string = subprocess.check_output([IgBlastCache.binary, '-version'])
        self._hasher = hashlib.md5(version_string)
        cache_home = os.environ.get('XDG_CACHE_HOME', os.path.expanduser('~/.cache'))
        self._cachedir = os.path.join(cache_home, 'igdiscover')
        logger.info('Caching IgBLAST results in %r', self._cachedir)
        self._lock = multiprocessing.Lock()

    def _path(self, digest):
        """Return path to cache file given a digest"""
        return os.path.join(self._cachedir, digest[:2], digest) + '.txt.gz'

    def _load(self, digest):
        try:
            with gzip.open(self._path(digest), 'rt') as f:
                return f.read()
        except FileNotFoundError:
            return None

    def _store(self, digest, data):
        path = self._path(digest)
        dir = os.path.dirname(path)
        os.makedirs(dir, exist_ok=True)
        with gzip.open(path, 'wt') as f:
            f.write(data)

    def retrieve(self, variable_arguments, fixed_arguments, blastdb_dir, fasta_str) -> str:
        hasher = self._hasher.copy()
        hasher.update(' '.join(fixed_arguments).encode())
        for gene in 'V', 'D', 'J':
            with open(os.path.join(blastdb_dir, gene + '.fasta'), 'rb') as f:
                hasher.update(f.read())
        hasher.update(fasta_str.encode())
        digest = hasher.hexdigest()
        data = self._load(digest)
        if data is None:
            with tempfile.TemporaryDirectory() as tmpdir:
                path = os.path.join(tmpdir, 'igblast.txt')
                full_arguments = [IgBlastCache.binary] + variable_arguments + fixed_arguments\
                    + ['-out', path]
                output = subprocess.check_output(full_arguments, input=fasta_str,
                    universal_newlines=True)
                assert output == ''
                with open(path) as f:
                    data = f.read()
            with self._lock:  # TODO does this help?
                self._store(digest, data)

        return data


def run_igblast(sequences, blastdb_dir, species, sequence_type, penalty=None, use_cache=True) -> str:
    """
    Run the igblastn command-line program.

    sequences -- list of Sequence objects
    blastdb_dir -- directory that contains BLAST databases. Files in that
    directory must be databases created by the makeblastdb program and have
    names V, D, and J.

    Return IgBLAST’s raw output as a string.
    """
    if sequence_type not in ('Ig', 'TCR'):
        raise ValueError('sequence_type must be "Ig" or "TCR"')
    variable_arguments = []
    for gene in 'V', 'D', 'J':
        variable_arguments += ['-germline_db_{gene}'.format(gene=gene),
            os.path.join(blastdb_dir, '{gene}'.format(gene=gene))]
    # An empty .aux suppresses a warning from IgBLAST. /dev/null does not work.
    empty_aux_path = pkg_resources.resource_filename('igdiscover', 'empty.aux')
    variable_arguments += ['-auxiliary_data', empty_aux_path]
    arguments = []
    if penalty is not None:
        arguments += ['-penalty', str(penalty)]
    if species is not None:
        arguments += ['-organism', species]

    arguments += [
        '-ig_seqtype', sequence_type,
        '-num_threads', '1',
        '-domain_system', 'imgt',
        '-num_alignments_V', '1',
        '-num_alignments_D', '1',
        '-num_alignments_J', '1',
        '-outfmt', '7 sseqid qstart qseq sstart sseq pident slen evalue',
        '-query', '-',
    ]
    fasta_str = ''.join(">{}\n{}\n".format(r.name, r.sequence) for r in sequences)

    if use_cache:
        global _igblastcache
        return _igblastcache.retrieve(variable_arguments, arguments, blastdb_dir, fasta_str)
    else:
        # For some reason, it has become unreliable to let IgBLAST 1.10 write its result
        # to standard output using "-out -". The data becomes corrupt. This does not occur
        # when calling igblastn in the shell and using the same syntax, only with
        # subprocess.check_output. As a workaround, we write the output to a temporary file.
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'igblast.txt')
            output = subprocess.check_output(['igblastn'] + variable_arguments + arguments
                + ['-out', path], input=fasta_str, universal_newlines=True)
            assert output == ''
            with open(path) as f:
                return f.read()


def chunked(iterable, chunksize: int):
    """
    Group the iterable into lists of length chunksize
    >>> list(chunked('ABCDEFG', 3))
    [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    """
    chunk = []
    for it in iterable:
        if len(chunk) == chunksize:
            yield chunk
            chunk = []
        chunk.append(it)
    if chunk:
        yield chunk


class Runner:
    """
    This is the target of a multiprocessing pool. The target needs to
    be pickleable, and because nested functions cannot be pickled,
    we need this separate class.

    It runs IgBLAST and parses the output for a list of sequences.
    """
    def __init__(self, blastdb_dir, species, sequence_type, penalty, database, use_cache):
        self.blastdb_dir = blastdb_dir
        self.species = species
        self.sequence_type = sequence_type
        self.penalty = penalty
        self.database = database
        self.use_cache = use_cache

    def __call__(self, sequences):
        """
        Return tuples (igblast_result, records) where igblast_result is the raw IgBLAST output
        and records is a list of (Extended-)IgBlastRecord objects (the parsed output).
        """
        igblast_result = run_igblast(sequences, self.blastdb_dir, self.species, self.sequence_type,
            self.penalty, self.use_cache)
        sio = StringIO(igblast_result)
        parser = IgBlastParser(sequences, sio, self.database)
        records = list(parser)
        assert len(records) == len(sequences)
        return igblast_result, records


def makeblastdb(fasta, database_name, prefix=''):
    """
    prefix -- prefix to add to sequence ids
    """
    n = 0
    with dnaio.open(fasta) as fr, open(database_name + ".fasta", "w") as db:
        for record in fr:
            name = prefix + record.name.split(maxsplit=1)[0]
            db.write(">{}\n{}\n".format(name, record.sequence))
            n += 1
    if n == 0:
        raise ValueError("FASTA file {} is empty".format(fasta))

    process_output = subprocess.check_output(
        ['makeblastdb', '-parse_seqids', '-dbtype', 'nucl', '-in', database_name + ".fasta", '-out',
            database_name],
        stderr=subprocess.STDOUT
    )
    if b'Error: ' in process_output:
        raise subprocess.SubprocessError()


class Database:
    def __init__(self, path, sequence_type):
        """path -- path to database directory with V.fasta, D.fasta, J.fasta"""
        self.path = path
        self.sequence_type = sequence_type
        self._v_records = self._read_fasta(os.path.join(path, 'V.fasta'))
        self.v = self._records_to_dict(self._v_records)
        self._j_records = self._read_fasta(os.path.join(path, 'J.fasta'))
        self.j = self._records_to_dict(self._j_records)
        self._cdr3_starts = dict()
        self._cdr3_ends = dict()
        for chain in ('heavy', 'kappa', 'lambda', 'alpha', 'beta', 'gamma', 'delta'):
            self._cdr3_starts[chain] = {name: cdr3_start(s, chain) for name, s in self.v.items()}
            self._cdr3_ends[chain] = {name: cdr3_end(s, chain) for name, s in self.j.items()}
        self.v_regions_nt, self.v_regions_aa = self._find_v_regions()

    @staticmethod
    def _read_fasta(path):
        records = []
        with dnaio.open(path) as sr:
            for record in sr:
                record.name = record.name.split(maxsplit=1)[0]
                records.append(record)
        return records

    @staticmethod
    def _records_to_dict(records):
        return {record.name: record.sequence.upper() for record in records}

    def v_cdr3_start(self, gene, chain):
        return self._cdr3_starts[chain][gene]

    def j_cdr3_end(self, gene, chain):
        return self._cdr3_ends[chain][gene]

    def _find_v_regions(self):
        """
        Run IgBLAST on the V sequences to determine the nucleotide and amino-acid sequences of the
        FR1, CDR1, FR2, CDR2 and FR3 regions
        """
        v_regions_nt = dict()
        v_regions_aa = dict()
        for record in igblast(self.path, self._v_records, self.sequence_type, threads=1):
            nt_regions = dict()
            aa_regions = dict()
            for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3'):
                nt_seq = record.region_sequence(region)
                if nt_seq is None:
                    break
                if len(nt_seq) % 3 != 0:
                    logger.warning('Length %s of %s region in %r is not divisible by three; region '
                        'info for this gene will not be available',
                        len(nt_seq), region, record.query_name)
                    # not codon-aligned, skip entire record
                    break
                nt_regions[region] = nt_seq
                try:
                    aa_seq = nt_to_aa(nt_seq)
                except ValueError as e:
                    logger.warning('The %s region could not be converted to amino acids: %s',
                        region, str(e))
                    break
                if '*' in aa_seq:
                    logger.warning('The %s region in %r contains a stop codon (%r); region info '
                        'for this gene will not be available',
                        region, record.query_name, aa_seq)
                    break
                aa_regions[region] = aa_seq
            else:
                v_regions_nt[record.query_name] = nt_regions
                v_regions_aa[record.query_name] = aa_regions

        return v_regions_nt, v_regions_aa


def igblast(database, sequences, sequence_type, species=None, threads=1, penalty=None,
        raw_output=None, use_cache=False):
    """
    Run IgBLAST, parse results and yield (Extended-)IgBlastRecord objects.

    database -- Path to database directory with V./D./J.fasta files *or* a Database object.
        If it is a path, then only IgBlastRecord objects are returned.
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    threads -- number of threads.
    raw_output -- If not None, raw IgBLAST output is written to this file
    """
    if isinstance(database, str):
        database_dir = database
        database = None
    else:
        database_dir = database.path
    with ExitStack() as stack:
        # Create the three BLAST databases in a temporary directory
        blastdb_dir = stack.enter_context(tempfile.TemporaryDirectory())
        for gene in ['V', 'D', 'J']:
            makeblastdb(
                os.path.join(database_dir, gene + '.fasta'),
                os.path.join(blastdb_dir, gene),
                prefix='%',
            )

        chunks = chunked(sequences, chunksize=1000)
        runner = Runner(blastdb_dir, species=species, sequence_type=sequence_type, penalty=penalty,
            database=database, use_cache=use_cache)
        pool = stack.enter_context(multiprocessing.Pool(threads) if threads > 1 else SerialPool())
        for igblast_output, igblast_records in pool.imap(runner, chunks, chunksize=1):
            if raw_output:
                raw_output.write(igblast_output)
            yield from igblast_records


def main(args):
    config = GlobalConfig()
    use_cache = config.use_cache
    if args.cache is not None:
        use_cache = args.cache
    if use_cache:
        global _igblastcache
        _igblastcache = IgBlastCache()
        logger.info('IgBLAST cache enabled')
    if args.threads == 0:
        args.threads = available_cpu_count()
    logger.info("Running IgBLAST on database sequences to find CDR/FR region locations")
    database = Database(args.database, args.sequence_type)
    logger.info("Running IgBLAST on input reads")
    detected_cdr3s = 0
    writer = TableWriter(sys.stdout)
    start_time = time.time()
    last_status_update = 0
    with ExitStack() as stack:
        if args.raw:
            raw_output = stack.enter_context(xopen(args.raw, 'w'))
        else:
            raw_output = None
        sequences = stack.enter_context(dnaio.open(args.fasta))
        sequences = islice(sequences, 0, args.limit)

        n = 0  # number of records processed so far
        for record in igblast(database, sequences, sequence_type=args.sequence_type,
                species=args.species, threads=args.threads, penalty=args.penalty,
                raw_output=raw_output, use_cache=use_cache):
            n += 1
            if args.rename is not None:
                record.query_name = "{}seq{}".format(args.rename, n)
            d = record.asdict()
            if d['CDR3_aa']:
                detected_cdr3s += 1
            try:
                writer.write(d)
            except IOError as e:
                if e.errno == errno.EPIPE:
                    sys.exit(1)
                raise
            if n % 1000 == 0:
                elapsed = time.time() - start_time
                if elapsed >= last_status_update + 60:
                    logger.info(
                        'Processed {:10,d} sequences at {:.3f} ms/sequence'.format(n, elapsed / n * 1E3))
                    last_status_update = elapsed
    elapsed = time.time() - start_time
    logger.info('Processed {:10,d} sequences at {:.1f} ms/sequence'.format(n, elapsed / n * 1E3))

    logger.info('%d IgBLAST assignments parsed and written', n)
    logger.info('CDR3s detected in %.1f%% of all sequences', detected_cdr3s / n * 100)
    if args.stats:
        stats = {'total': n, 'detected_cdr3s': detected_cdr3s}
        with open(args.stats, 'w') as f:
            json.dump(stats, f)
            print(file=f)
