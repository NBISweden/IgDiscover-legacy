"""
Merge paired-end reads (a wrapper around PEAR)

This script can also manage a cache of already-merged files.
"""
import os
import sys
from pathlib import Path
import logging
import tempfile
import subprocess
import hashlib
import shutil

from ..config import GlobalConfig

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--threads', '-j', default=None, type=int, help='Number of threads')
    arg('--no-cache', default=None, action='store_true',
        help='Disable cache. Default: Determined by configuration')
    arg('reads1', help='Forward reads FASTQ file')
    arg('reads2', help='Reverse reads FASTQ file')
    arg('output',  help='Output file (compressed FASTQ)')


def compute_hash(path1: str, path2: str):
    result = subprocess.run(['pear', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pear_version = None
    for line in result.stdout.split(b'\n'):
        if line.startswith(b'PEAR'):
            pear_version = line
            break
    assert pear_version is not None

    hasher = hashlib.md5(pear_version)
    for path in path1, path2:
        with open(path, 'rb') as f:
            while True:
                chunk = f.read(1048576)
                if not chunk:
                    break
                hasher.update(chunk)
    return hasher.hexdigest()


def run_pear(path1: str, path2: str, output: str, log_output=None, threads: int=None):
    gzip = 'pigz' if shutil.which('pigz') else 'gzip'
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        args = ['pear']
        if threads is not None:
            args += ['-j', str(threads)]
        args += [
            '-f', path1,
            '-r', path2,
            '-o', tmpdir / 'merged'
        ]
        if log_output is None:
            subprocess.run(args, check=True)
        else:
            result = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with open(log_output, 'wb') as f:
                f.write(result.stdout)

        subprocess.run([gzip, tmpdir / 'merged.assembled.fastq'], check=True)
        shutil.move(tmpdir / 'merged.assembled.fastq.gz', output)


def run_pear_cached(path1: str, path2: str, output: str, threads: int = None):
    cache_home = Path(os.environ.get('XDG_CACHE_HOME', os.path.expanduser('~/.cache')))
    cachedir = cache_home / 'igdiscover' / 'mergedreads'
    cachedir.mkdir(parents=True, exist_ok=True)

    md5sum = compute_hash(path1, path2)
    cached_merged = cachedir / (md5sum + '.fastq.gz')
    cached_log = cachedir / (md5sum + '.log')
    if not cached_merged.exists() or not cached_log.exists():
        run_pear(path1, path2, output=cached_merged, log_output=cached_log, threads=threads)
        logger.info('PEAR result copied to cache')
    else:
        logger.info('PEAR result found in cache:')
        logger.info(cached_merged)
        logger.info(cached_log)
    shutil.copy(cached_merged, output)
    with open(cached_log) as f:
        sys.stdout.write(f.read())

    # Update mtime of the cache file. Let’s us find unused cache entries.
    for path in cached_merged, cached_log:
        try:
            path.touch()
        except PermissionError:
            pass


def main(args):
    use_cache = GlobalConfig().use_cache
    if args.no_cache:
        use_cache = False
    if use_cache:
        logger.info('Cache enabled')
        func = run_pear_cached
    else:
        func = run_pear
    func(path1=args.reads1, path2=args.reads2, output=args.output, threads=args.threads)
