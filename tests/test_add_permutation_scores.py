"""test_add_permutation_scores.py - exact and drawn p, tables and replicates

input:   small ACGT tables, one to enumerate and two deep enough to draw
output:  pytest assertions
usage:   pytest tests/test_add_permutation_scores.py
env:     numpy and scipy, reached through the module under test
notes:   the pmf takes integer counts while the scorer takes floats, and the
         enumerated path is where the two meet
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from add_permutation_scores import (
    MIN_HIGHER_SCORING,
    REDRAW,
    counts_for,
    enumerate_splits,
    enumerate_tables,
    fit,
    permutation_test,
    permute_reads,
    pseudocounts,
    redraw,
    scorer_for,
    single_replicate_loglik,
    table_pmf,
)

COUNTS1 = np.array([[8, 1, 0, 1], [30, 2, 1, 1], [400, 60, 20, 20]],
                   dtype=float)[:, None]
COUNTS2 = np.array([[7, 2, 1, 0], [28, 9, 1, 2], [380, 20, 30, 70]],
                   dtype=float)[:, None]
DRAWS = 50
SCORER = scorer_for("mismatch", 20)


@pytest.fixture(scope="module")
def result():
    return permutation_test((COUNTS1, COUNTS2), SCORER, seed=42, draws=DRAWS)


def test_p_is_exact_up_to_draws_possible_tables(result):
    assert result["p_exact"].tolist() == [True, False, False]


def test_every_p_is_a_probability(result):
    assert np.all(np.isfinite(result["pvalue"]))
    assert np.all((result["pvalue"] > 0) & (result["pvalue"] <= 1))


def test_the_margin_is_the_score_above_the_null(result):
    assert result["score_minus_max_null"][2] > 0
    assert np.all(result["score_minus_max_null"] >= 0)


def test_a_site_draws_until_enough_score_higher_or_the_maximum():
    # permutations at or above the observed score: 0 of 5000 at site 0, and
    # 2 and 13 of the first 50 at sites 1 and 2
    sites = (np.array([[[400, 0, 0, 0]], [[100, 10, 5, 5]], [[100, 10, 5, 5]]],
                      dtype=float),
             np.array([[[0, 0, 0, 400]], [[85, 25, 5, 5]], [[90, 20, 5, 5]]],
                      dtype=float))
    results = {"mismatch": permutation_test(sites, SCORER, seed=42,
                                            draws=DRAWS)}
    first = results["mismatch"]["pvalue"][2]
    maximum = DRAWS * REDRAW ** 2
    args = SimpleNamespace(max_permutations=maximum, seed=42, qual=20,
                           min_higher_scoring=MIN_HIGHER_SCORING)
    redraw(results, {"mismatch": sites}, args, map)
    n, p = (results["mismatch"][k] for k in ("n_permutations", "pvalue"))
    assert n[0] == maximum and p[0] == 1 / (maximum + 1)
    assert DRAWS < n[1] < maximum and p[1] == MIN_HIGHER_SCORING / n[1]
    assert n[2] == DRAWS and p[2] == first


def test_a_site_with_permutations_scoring_higher_keeps_the_first_draws():
    same = np.array([[[500, 50, 50, 50]]], dtype=float)
    found = permutation_test((same, same), SCORER, seed=42, draws=DRAWS)
    assert found["n_permutations"][0] == DRAWS
    n_higher = found["pvalue"][0] * (DRAWS + 1) - 1
    assert n_higher > 0 and np.isclose(n_higher, round(n_higher))


def test_replicates_stay_separate_through_the_test():
    counts1 = np.array([[[400, 60, 20, 20], [390, 50, 30, 30]]], dtype=float)
    counts2 = np.array([[[380, 20, 30, 70], [370, 25, 35, 70]]], dtype=float)
    found = permutation_test((counts1, counts2), SCORER, seed=42, draws=DRAWS)
    assert 0 < found["pvalue"][0] <= 1


def test_a_permutation_keeps_every_depth_and_the_site_counts():
    counts = np.array([[[40, 6, 2, 2], [30, 5, 3, 2], [20, 9, 1, 0]]])
    permuted = permute_reads(counts, np.random.default_rng(1), 200)
    assert (permuted.sum(-1) == counts.sum(-1)[:, None]).all()
    assert (permuted.sum(2) == counts.sum(1)[:, None]).all()


def test_two_group_tables_are_the_splits_and_pmf_sums_to_one():
    pooled = np.array([15, 3, 1, 1])
    two = enumerate_tables(pooled, [10, 10], 10_000)
    np.testing.assert_array_equal(two[:, 0], enumerate_splits(pooled, 10,
                                                              10_000))
    three = enumerate_tables(pooled, [8, 7, 5], 10_000)
    assert np.isclose(table_pmf(three, pooled).sum(), 1)


def test_single_replicate_loglik_is_the_limit_the_stepwise_fit_approaches():
    counts = np.random.default_rng(2).integers(0, 300, size=(200, 1, 4))
    data = pseudocounts(counts.astype(float), 20)
    limit = single_replicate_loglik(data)[0]
    gap = limit - fit(data)[1]
    assert np.all(gap >= 0) and np.all(gap < 2.1e-3)


@pytest.mark.parametrize("error_type", ["mismatch", "deletion", "insertion"])
def test_each_replicate_is_one_row(error_type):
    table = pd.DataFrame({
        "bases11": ["8,1,0,1"], "bases12": ["2,0,1,0"], "bases21": ["7,2,1,0"],
        "deletions11": [3], "deletions12": [1], "deletions21": [2],
        "insertions11": [1], "insertions12": [2], "insertions21": [0]})
    first, second = counts_for(table, error_type)
    width = 4 if error_type == "mismatch" else 2
    assert first.shape == (1, 2, width) and second.shape == (1, 1, width)
    if error_type == "mismatch":
        np.testing.assert_array_equal(first[0, 1], [2, 0, 1, 0])
    else:
        found = 1 if error_type == "deletion" else 2
        assert first[0, 1].tolist() == [found, 3 + 1 - found]
