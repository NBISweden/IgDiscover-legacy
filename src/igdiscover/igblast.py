"""
This provides functions for running the "igblastn" command-line tool
"""
import csv
import os
import shlex
import multiprocessing
import subprocess
from contextlib import ExitStack
from dataclasses import dataclass
from io import StringIO
import hashlib
from typing import Dict

import pkg_resources
import logging
import tempfile
import gzip

import dnaio

from .utils import SerialPool, nt_to_aa
from .species import cdr3_start, cdr3_end


logger = logging.getLogger(__name__)


def escape_shell_command(command):
    return " ".join(shlex.quote(arg) for arg in command)


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
        '-outfmt', '19',  # AIRR format
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
            command = ['igblastn'] + variable_arguments + arguments + ['-out', path]
            logger.debug("Running %s", escape_shell_command(command))
            output = subprocess.check_output(command, input=fasta_str, universal_newlines=True)
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


class RawRunner:
    """
    This is the target of a multiprocessing pool. The target needs to
    be pickleable, and because nested functions cannot be pickled,
    we need this separate class.

    It runs IgBLAST and returns raw AIRR-formatted output
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
        )


class MakeBlastDbError(subprocess.CalledProcessError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __str__(self):
        return f"Running '{escape_shell_command(self.cmd)}' failed with " \
               f"exit code {self.returncode}. " \
               f"Standard output:\n{self.output.decode()}\n" \
               f"Standard error:\n{self.stderr.decode()}"


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

    command = [
        'makeblastdb', '-parse_seqids', '-dbtype', 'nucl', '-in', database_name + ".fasta", '-out',
        database_name
    ]
    logger.debug("Running %s", escape_shell_command(command))
    try:
        process_output = subprocess.check_output(command, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise MakeBlastDbError(e.returncode, command, output=e.output, stderr=e.stderr) from None

    if b'Error: ' in process_output:
        raise MakeBlastDbError(0, command, stderr=process_output) from None


class Database:
    def __init__(self, path, sequence_type):
        """path -- path to database directory with V.fasta, D.fasta, J.fasta"""
        self.path = path
        self.sequence_type = sequence_type
        self._v_records = self._read_fasta(os.path.join(path, 'V.fasta'))
        self.v = self._records_to_dict(self._v_records)
        self._d_records = self._read_fasta(os.path.join(path, 'D.fasta'))
        self.d = self._records_to_dict(self._d_records)
        self._j_records = self._read_fasta(os.path.join(path, 'J.fasta'))
        self.j = self._records_to_dict(self._j_records)
        self._cdr3_starts = dict()
        self._cdr3_ends = dict()
        for locus in ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"):
            self._cdr3_starts[locus] = {name: cdr3_start(s, locus) for name, s in self.v.items()}
            self._cdr3_ends[locus] = {name: cdr3_end(s, locus) for name, s in self.j.items()}
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

    def v_cdr3_start(self, gene, locus):
        return self._cdr3_starts[locus][gene]

    def j_cdr3_end(self, gene, locus):
        return self._cdr3_ends[locus][gene]

    def _find_v_regions(self):
        """
        Run IgBLAST on the V sequences to determine the nucleotide and amino-acid sequences of the
        FR1, CDR1, FR2, CDR2 and FR3 regions
        """
        v_regions_nt = dict()
        v_regions_aa = dict()
        for record in igblast_records(self.path, self._v_records, self.sequence_type):
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


def igblast_records(database, sequences, sequence_type, species=None, penalty=None,
        use_cache=False):
    """
    Run IgBLAST, parse results and yield RegionRecords objects.

    database -- Path to database directory with V./D./J.fasta files
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    """
    # Create the BLAST databases in a temporary directory
    with tempfile.TemporaryDirectory() as blastdb_dir:
        make_vdj_blastdb(blastdb_dir, database)
        igblast_result = run_igblast(sequences, blastdb_dir, species, sequence_type, penalty, use_cache)
        sio = StringIO(igblast_result)
        for record in parse_region_records(sio):
            yield record


def parse_region_records(file):
    csv.register_dialect(
        "airr",
        delimiter="\t",
        lineterminator="\n",
        strict=True,
    )
    reader = csv.DictReader(file, dialect="airr")
    for record in reader:
        yield RegionsRecord(
            query_name=record["sequence_id"],
            fields=record,
        )


@dataclass
class RegionsRecord:
    query_name: str
    fields: Dict[str, str]

    COLUMNS_MAP = {
        "CDR1": "cdr1",
        "CDR2": "cdr2",
        "FR1": "fwr1",
        "FR2": "fwr2",
        "FR3": "fwr3",
    }

    def region_sequence(self, region: str) -> str:
        """
        Return the nucleotide sequence of a named region. Allowed names are:
        CDR1, CDR2, CDR3, FR1, FR2, FR3. Sequences are extracted from the full read
        using begin and end coordinates from IgBLAST’s "alignment summary" table.
        """
        if region not in self.COLUMNS_MAP:
            raise KeyError(f"Region '{region}' not allowed")
        return self.fields[self.COLUMNS_MAP[region]]


def igblast_parallel_chunked(
    database,
    sequences,
    sequence_type,
    species=None,
    threads=1,
    penalty=None,
    use_cache=False,
):
    """
    Run IgBLAST on the input sequences and yield AIRR-formatted results.

    The input is split up into chunks of 1000 sequences and distributed to
    *threads* number of IgBLAST instances that run in parallel.

    database -- Path to database directory with V./D./J.fasta files
    sequences -- an iterable of Sequence objects
    sequence_type -- 'Ig' or 'TCR'
    threads -- number of threads.
    """
    with ExitStack() as stack:
        # Create the three BLAST databases in a temporary directory
        blastdb_dir = stack.enter_context(tempfile.TemporaryDirectory())
        make_vdj_blastdb(blastdb_dir, database)

        chunks = chunked(sequences, chunksize=1000)
        runner = RawRunner(
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


def make_vdj_blastdb(blastdb_dir, database_dir):
    """Run makeblastdb for all {V,D,J}.fasta in the database_dir"""

    for gene in ["V", "D", "J"]:
        # Without adding the "%" prefix, IgBLAST reports record names that look like GenBank
        # ids as "gb|original_name|".
        makeblastdb(
            os.path.join(database_dir, gene + ".fasta"),
            os.path.join(blastdb_dir, gene),
            prefix="%",
        )
