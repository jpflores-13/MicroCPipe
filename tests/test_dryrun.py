"""Snakemake dry runs of both workflows against the tiny fixture in tests/dryrun/.

A dry run parses the samplesheet and config, builds the DAG and checks that
input files exist; it runs no tools. Skipped when snakemake is not installed.
Run from anywhere: the tests chdir to the repository root.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFIG = "tests/dryrun/config.yaml"

pytestmark = pytest.mark.skipif(shutil.which("snakemake") is None, reason="snakemake not installed")


def dry_run(snakefile, *extra):
    cmd = ["snakemake", "-n", "-s", f"workflows/{snakefile}", "--cores", "1",
           "--configfile", "config/MicroCConfig.yaml", CONFIG, *extra]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)


def test_core_dry_run():
    r = dry_run("MicroCPipe.snakefile")
    assert r.returncode == 0, r.stderr
    # three samplesheet rows -> three alignment jobs; two groups (Tech_Rep 1 and 2)
    assert "align" in r.stdout and "dedup" in r.stdout and "qc_report" in r.stdout
    assert "CI_CellA_1_1" in r.stdout and "CI_CellA_1_2" in r.stdout


def test_core_counts():
    r = dry_run("MicroCPipe.snakefile")
    rows = {ln.split()[0]: int(ln.split()[1]) for ln in r.stdout.splitlines()
            if ln.split() and ln.split()[0] in {"align", "dedup", "cool", "mcool", "hic"}
            and len(ln.split()) > 1 and ln.split()[1].isdigit()}
    assert rows.get("align") == 3
    assert rows.get("dedup") == 2


def test_merge_dry_run():
    r = dry_run("MicroCMerge.snakefile")
    assert r.returncode == 0, r.stderr
    assert "merge_pairs" in r.stdout and "merged_CI_CellA_1" in r.stdout
    assert "merged_CI_CellA_1.mcool" in r.stdout


def test_bad_output_format_is_rejected():
    r = dry_run("MicroCPipe.snakefile", "--config", "outputs=['bam']")
    assert r.returncode != 0
    assert "unknown format" in (r.stderr + r.stdout)
