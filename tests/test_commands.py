import os
import sys
import contextlib
import shutil
from pathlib import Path

from xopen import xopen
import pytest

from igdiscover.__main__ import main
from .utils import datapath, resultpath, files_equal, convert_fastq_to_fasta
from igdiscover.cli.init import run_init
from igdiscover.cli.config import print_configuration, modify_configuration
from igdiscover.cli.run import run_snakemake
from igdiscover.cli.clonotypes import run_clonotypes


@pytest.fixture
def run(tmpdir):
    def _run(args, expected):
        """
        Run IgDiscover, redirecting stdout to a temporary file.
        Then compare the output with the contents of an expected file.
        """
        outpath = str(tmpdir.join('output'))
        print('Running:', ' '.join(args))
        with open(outpath, 'w') as f:
            old_stdout = sys.stdout
            sys.stdout = f
            main(args)
            sys.stdout = old_stdout
        assert files_equal(expected, outpath)

    return _run


@pytest.fixture
def pipeline_dir(tmp_path):
    """An initialized pipeline directory"""
    pipeline_path = tmp_path / "initializedpipeline"
    init_testdata(pipeline_path)
    return pipeline_path


def init_testdata(directory):
    run_init(
        database="testdata/database",
        reads1="testdata/reads.1.fastq.gz",
        directory=str(directory),
    )
    with chdir(directory):
        modify_configuration([("barcode_length_3prime", "21")])


@contextlib.contextmanager
def chdir(path):
    previous_path = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(previous_path)


@pytest.fixture(scope="session")
def filtered_tab_session(tmp_path_factory):
    """Generate iteration-01/filtered.tsv.gz"""

    pipeline_dir = tmp_path_factory.mktemp("pipedir") / "pipedir"
    init_testdata(pipeline_dir)
    with chdir(pipeline_dir):
        run_snakemake(targets=["iteration-01/filtered.tsv.gz"])
    return pipeline_dir


@pytest.fixture
def has_filtered_tsv(filtered_tab_session, tmp_path):
    """
    Give a fresh copy of a pipeline dir in which iteration-01/filtered.tsv.gz
    is guaranteed to exist
    """
    pipeline_dir = tmp_path / "has_filtered_tsv"
    shutil.copytree(
        filtered_tab_session,
        pipeline_dir,
        symlinks=True,
        ignore=shutil.ignore_patterns((".snakemake")),
    )
    return pipeline_dir


def test_main():
    with pytest.raises(SystemExit) as exc:
        main(['--version'])
    assert exc.value.code == 0


def test_group_by_barcode_only(run):
    args = ['group', '-b', '4', datapath('ungrouped.fasta')]
    run(args,  resultpath('grouped-by-barcode-only.fasta'))


def test_group_by_pseudo_cdr3(run):
    args = ['group', '-b', '4', '--pseudo-cdr3=-5:-2', '--trim-g', datapath('ungrouped.fasta')]
    run(args,  resultpath('grouped.fasta'))


def test_group_by_pseudo_cdr3_barcode_at_end(run):
    args = ['group', '-b', '-4', '--pseudo-cdr3=1:3', datapath('ungrouped.fasta')]
    run(args, resultpath('grouped2.fasta'))


def test_clusterplot(tmpdir):
    main(['clusterplot', '-m', '10', datapath('clusterplot.tab.gz'), str(tmpdir)])
    assert tmpdir.join('IGHV1-1801.png').check()


def test_igblastwrap(run):
    args = ['igblastwrap', '--threads=1', datapath('database/'), datapath('igblast.fasta')]
    run(args, resultpath('assigned.tsv'))


def test_run_init(pipeline_dir):
    assert pipeline_dir.is_dir()
    assert (pipeline_dir / "igdiscover.yaml").exists()


def test_print_configuration(pipeline_dir):
    print_configuration(path=pipeline_dir / "igdiscover.yaml")


def test_modify_configuration(pipeline_dir):
    modify_configuration(
        settings=[("d_coverage", "12"), ("j_discovery.allele_ratio", "0.37")],
        path=str(pipeline_dir / "igdiscover.yaml"),
    )
    from ruamel.yaml import YAML
    yaml = YAML(typ="safe", pure=True)
    with open(pipeline_dir / "igdiscover.yaml") as f:
        config = yaml.load(f)
    assert config["d_coverage"] == 12
    assert config["j_discovery"]["allele_ratio"] == 0.37


def test_dryrun(pipeline_dir):
    with chdir(pipeline_dir):
        run_snakemake(dryrun=True)


def test_primers(pipeline_dir):
    # Test whether specifying primer sequences leads to a SyntaxError
    with chdir(pipeline_dir):
        modify_configuration(
            settings=[
                ("forward_primers", "['CGTGA']"),
                ("reverse_primers", "['TTCAC']"),
            ],
        )
        run_snakemake(dryrun=True)


# TODO slow
def test_only_forward_primer(pipeline_dir):
    # issue #107 (broken symlink)
    with chdir(pipeline_dir):
        modify_configuration(settings=[("forward_primers", "['CGTGA']")])
        # Create some dummy files so we don’t need to run irrelevant steps of the pipeline
        r = Path("reads")
        r.mkdir()
        with xopen(r / "2-merged.fastq.gz", "w") as f:
            pass
        s = Path("stats")
        s.mkdir()
        with open(s / "merging-successful", "w") as f:
            pass
        with open(s / "reads.json", "w") as f:
            f.write('{"total": 0}')
        with open(s / "trimmed.json", "w") as f:
            f.write('{"trimmed": 0}')
        run_snakemake(targets=["reads/sequences.fasta.gz"])


def test_flash(pipeline_dir):
    # Test using FLASH and parsing its log output
    with chdir(pipeline_dir):
        modify_configuration(settings=[("merge_program", "flash")])
        run_snakemake(targets=["stats/reads.json"])
        # Ensure FLASH was actually run
        assert (pipeline_dir / "reads/2-flash.log").exists()


def test_snakemake_assigned_tab(has_filtered_tsv):
    assert (has_filtered_tsv / "iteration-01/filtered.tsv.gz").exists()
    assert not (has_filtered_tsv / "iteration-01/new_V_germline.tab").exists()


def test_snakemake_exact_tab(has_filtered_tsv):
    with chdir(has_filtered_tsv):
        run_snakemake(targets=["iteration-01/exact.tab"])
    assert (has_filtered_tsv / "iteration-01/exact.tab").exists()


def test_snakemake_final(has_filtered_tsv):
    with chdir(has_filtered_tsv):
        run_snakemake(targets=["nofinal"])
    assert (has_filtered_tsv / "iteration-01/new_V_germline.tab").exists()
    assert not (has_filtered_tsv / "final/assigned.tsv.gz").exists()

    with chdir(has_filtered_tsv):
        run_snakemake()
    assert (has_filtered_tsv / "final/assigned.tsv.gz").exists()


def test_clonotypes(has_filtered_tsv):
    run_clonotypes(has_filtered_tsv / "iteration-01/assigned.tsv.gz", limit=5)


def test_fastq_input(has_filtered_tsv, tmp_path):
    # Use merged reads from already-run pipeline as input for a new run
    single_reads = has_filtered_tsv / "reads" / "2-merged.fastq.gz"
    directory = tmp_path / "singleend-fastq"
    run_init(
        database="testdata/database",
        single_reads=str(single_reads),
        directory=str(directory),
    )
    with chdir(directory):
        modify_configuration([("barcode_length_3prime", "21")])
        run_snakemake(targets=["stats/reads.json"])


def test_fasta_input(has_filtered_tsv, tmp_path):
    fasta_path = tmp_path / "justfasta.fasta"
    convert_fastq_to_fasta(
        has_filtered_tsv / "reads" / "2-merged.fastq.gz",
        fasta_path,
    )
    directory = tmp_path / "singleend-fasta"
    run_init(
        database="testdata/database",
        single_reads=str(fasta_path),
        directory=str(directory),
    )
    with chdir(directory):
        modify_configuration([("barcode_length_3prime", "21")])
        run_snakemake(targets=["stats/reads.json"])
