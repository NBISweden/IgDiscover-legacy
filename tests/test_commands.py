import sys
import pytest

from igdiscover.__main__ import main
from .utils import datapath, resultpath, files_equal
from igdiscover.cli.init import run_init
from igdiscover.cli.config import print_configuration, modify_configuration


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
    pipelinedir = tmp_path / "testdir"
    run_init(
        database="testdata/database",
        reads1="testdata/reads.1.fastq.gz",
        directory=str(pipelinedir),
    )
    return pipelinedir


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


def test_igblast(run):
    args = ['igblast', '--threads=1', datapath('database/'), datapath('igblast.fasta')]
    run(args, resultpath('assigned.tab'))


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
    import ruamel.yaml
    with open(pipeline_dir / "igdiscover.yaml") as f:
        config = ruamel.yaml.safe_load(f)
    assert config["d_coverage"] == 12
    assert config["j_discovery"]["allele_ratio"] == 0.37
