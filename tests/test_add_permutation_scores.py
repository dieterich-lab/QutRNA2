"""test_add_permutation_scores.py - both p value paths on three positions

input:   ACGT counts small enough to enumerate and one deep enough to draw
output:  pytest assertions
usage:   pytest tests/test_add_permutation_scores.py
env:     numpy and scipy, reached through the module under test
notes:   the pmf takes integer counts while the scorer takes floats, and the
         enumerated path is where the two meet
"""
import numpy as np
import pytest
from add_permutation_scores import permutation_test

COUNTS1 = np.array([[8, 1, 0, 1], [30, 2, 1, 1], [400, 60, 20, 20]],
                   dtype=float)
COUNTS2 = np.array([[7, 2, 1, 0], [28, 9, 1, 2], [380, 20, 30, 70]],
                   dtype=float)
DRAWS = 50


@pytest.fixture(scope="module")
def result():
    return permutation_test((COUNTS1, COUNTS2), qual=20, seed=42, draws=DRAWS)


def test_both_paths_run(result):
    assert result["p_exact"].tolist() == [True, True, False]


def test_every_p_is_a_probability(result):
    assert np.all(np.isfinite(result["pvalue"]))
    assert np.all((result["pvalue"] > 0) & (result["pvalue"] <= 1))


def test_a_drawn_p_stops_at_the_sampling_limit(result):
    assert result["pvalue"][2] >= 1 / (DRAWS + 1)


def test_the_margin_is_the_score_above_the_null(result):
    assert result["score_minus_max_null"][2] > 0
    assert np.all(result["score_minus_max_null"] >= 0)
