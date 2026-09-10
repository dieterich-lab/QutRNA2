"""Unit tests for workflow/scripts/add_scores.R.

JACUSA2 writes a header and no rows when no position reaches min_cov in every
replicate of both conditions. The error has to name that cause, since the
reader otherwise reports only a line count.
"""
import os
import shutil
import subprocess

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts",
                      "add_scores.R")

HEADER = ("#contig\tstart\tend\tname\tscore\tstrand\tbases11\tbases21\t"
          "filter\tref\n")


def test_header_without_rows_names_min_cov(tmp_path):
    assert shutil.which("Rscript"), \
        "Rscript is missing; run the tests inside tests/conda_env_tests.yaml"

    j2 = tmp_path / "JACUSA2.out"
    j2.write_text("## JACUSA2 Version: 2.1.16 (main) call-2\n" + HEADER)
    out = tmp_path / "scores.tsv"

    r = subprocess.run(
        ["Rscript", SCRIPT, "-s", "M::mismatch_score", "-o", str(out),
         str(j2)],
        capture_output=True, text=True)

    assert r.returncode != 0
    assert "no data rows" in r.stderr
    assert "min_cov" in r.stderr
    assert not out.exists()
