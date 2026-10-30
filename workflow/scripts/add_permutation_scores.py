"""add_permutation_scores.py - permutation p and BH FDR over a JACUSA2 call-2 run

input:   a JACUSA2 call-2 output file, two conditions, any replicates each
output:  a TSV of the JACUSA2 columns bar the subsampled ones, plus the
         permutation columns
usage:   python add_permutation_scores.py JACUSA2.out -o perm.tsv --seed 42
env:     python3 with numpy, scipy, pandas. JACUSA2 is not called
notes:   new columns are perm_<type>_<name>. type is mismatch, deletion or
         insertion. name is score, score_minus_max_null, pvalue, fdr,
         p_exact, n_splits or n_permutations. n_splits counts the possible
         read splits where p is exact, and is 0 where p is drawn
"""
import itertools
import os
import re
import sys
import time
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from contextlib import nullcontext
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.special import digamma, gammaln, polygamma
from scipy.stats import false_discovery_control, multivariate_hypergeom

CATEGORIES = 4
ESTIMATED_ERROR = 0.01
MAX_Q = 40
MAX_ITERATIONS = 100
EPSILON = 0.001
BACKTRACK_LAMBDA = 0.9
INDEL_PSEUDOCOUNT = 1.0

# one fit call over many positions pays the setup once, and bounds memory
BATCH_TABLES = 250_000
# a split scoring the same as the observed data counts as at or above it,
# whatever the float rounding
TIE_TOLERANCE = 1e-9
# positions per task. 25 gives even a short input more tasks than processes.
# not a flag, since each task's seed is the seed plus its task index
CHUNK = 25
# also the possible read splits up to which p is exact, and the permutations
# whose highest score score_minus_max_null subtracts
MIN_PERMUTATIONS = 1000
# a site draws until this many permutations score at or above its observed
# score, or until MAX_PERMUTATIONS. Besag and Clifford 1991
MIN_HIGHER_SCORING = 10
MAX_PERMUTATIONS = 100_000
SEED = 42
# each round multiplies a site's permutation count by this
REDRAW = 10

# the prefix of JACUSA2's count columns
ERROR_TYPES = {
    "mismatch": "bases",
    "deletion": "deletions",
    "insertion": "insertions",
}
ADDED_COLUMNS = ["score", "score_minus_max_null", "pvalue", "fdr",
                 "p_exact", "n_splits", "n_permutations"]


class InputError(Exception):
    """The input file cannot be scored, with the reason in the message."""


# ---------------------------------------------------------------- the score
# JACUSA2's Dirichlet-multinomial call-2 score, following
# lib/estimate/MinkaEstimateDirMultAlpha.java over all positions at once


# JACUSA2's Java, which this port is checked against, has the same cap
def error_prob(qual):
    """Convert a Phred score to an error probability. Caps at MAX_Q."""
    return 10.0 ** (-min(qual, MAX_Q) / 10.0)


def pseudocounts(counts, qual):
    """Turn raw base counts into the pseudo-counts the fit sees.

    A little error mass is spread from each observed base to the other three.
    Without it an unobserved base is zero, where logGamma(0) is undefined.
    """
    total = counts.sum(-1, keepdims=True)
    others = total - counts
    mean_error = np.divide(error_prob(qual) / (CATEGORIES - 1) * others, total,
                           out=np.zeros_like(counts), where=total > 0)

    return counts + (mean_error + ESTIMATED_ERROR) * others / (CATEGORIES - 1)


def loglik(alpha, data):
    """DirMult log-likelihood per position, summed over replicate rows.

    alpha is (N, K). data is (N, R, K). The result is (N,).
    """
    alpha_sum = alpha.sum(-1)[:, None]
    outer = gammaln(alpha_sum) - gammaln(data.sum(-1) + alpha_sum)
    inner = gammaln(data + alpha[:, None]) - gammaln(alpha)[:, None]

    return outer.sum(-1) + inner.sum((-1, -2))


def init_alpha(data, *, fallback):
    """Pick the starting alpha for a fit.

    The default is each category's pooled proportion, JACUSA2's MeanAlphaInit.
    The fallback is its smallest proportion across replicate rows, MinAlphaInit,
    which JACUSA2 reaches once the default has failed.
    """
    if fallback:
        return (data / data.sum(-1, keepdims=True)).min(1)

    summed = data.sum(1)

    return summed / summed.sum(-1, keepdims=True)


def _invalid(alpha):
    return (alpha < 0).any(-1) | ~np.isfinite(alpha).all(-1)


def _newton_step(alpha, data, *, backtrack):
    """Run one Minka iteration.

    Returns the proposed alpha, and a mask marking the positions that lack a
    valid alpha.
    """
    replicates = data.shape[1]
    alpha_sum = alpha.sum(-1)[:, None]
    row_sums = data.sum(-1)

    gradient = (replicates * digamma(alpha_sum)
                - digamma(row_sums + alpha_sum).sum(-1)[:, None]
                + (digamma(data + alpha[:, None])
                   - digamma(alpha)[:, None]).sum(1))
    q = (polygamma(1, data + alpha[:, None])
         - polygamma(1, alpha)[:, None]).sum(1)

    q_inverse_sum = (1.0 / q).sum(-1)
    z = (replicates * polygamma(1, alpha_sum[:, 0])
         - polygamma(1, row_sums + alpha_sum).sum(-1))
    b = (gradient / q).sum(-1) / (1.0 / z + q_inverse_sum)

    proposed = alpha - (gradient - b[:, None]) / q
    bad = _invalid(proposed)
    if backtrack and bad.any():
        damped = alpha[bad] - BACKTRACK_LAMBDA * (
            gradient[bad] - q_inverse_sum[bad][:, None]) / q[bad]
        proposed[bad] = damped
        bad[bad] = _invalid(damped)

    return proposed, bad


def fit(data, *, fallback=False):
    """Maximise the DirMult log-likelihood for every position at once.

    data is (N, R, K). Returns the fitted alpha, its log-likelihood, and a mask
    of the positions whose every step stayed valid.
    """
    alpha = init_alpha(data, fallback=fallback)
    active = np.ones(len(data), dtype=bool)
    converged = active.copy()

    for _ in range(MAX_ITERATIONS):
        idx = np.flatnonzero(active)
        if not idx.size:
            break

        proposed, bad = _newton_step(alpha[idx], data[idx], backtrack=fallback)
        ll_old = loglik(alpha[idx], data[idx])

        converged[idx[bad]] = False
        active[idx[bad]] = False

        good = idx[~bad]
        alpha[good] = proposed[~bad]
        delta = np.abs(loglik(alpha[good], data[good]) - ll_old[~bad])
        active[good[delta <= EPSILON]] = False

    return alpha, loglik(alpha, data), converged


def single_replicate_loglik(data):
    """The fit of a single replicate, (N, 1, K), in closed form.

    One count vector shows no overdispersion: the likelihood rises without
    bound in the precision, towards sum of x log(x / N). Returns that limit
    and an all-true convergence mask.
    """
    counts = data[:, 0]
    loglik = (counts * np.log(counts / counts.sum(-1, keepdims=True))).sum(-1)

    return loglik, np.ones(len(counts), dtype=bool)


def _with_replicates(counts):
    """Add a replicate axis where one is absent. Takes (N, K) or (N, R, K)."""
    counts = np.asarray(counts, dtype=float)

    return counts[:, None] if counts.ndim == 2 else counts


def log_likelihood_ratio(table1, table2):
    """Fit each condition and the pool, and return the likelihood ratio.

    Tables are (N, R, K) pseudo-counts. Returns the statistic and a mask of the
    positions where all three fits converged on the default init.
    """
    data = [table1, table2]
    data.append(np.concatenate(data, axis=1))

    logliks, fit_converged = zip(*(single_replicate_loglik(d) if d.shape[1] == 1
                                   else fit(d)[1:] for d in data))
    converged = np.logical_and.reduce(fit_converged)

    # JACUSA2 refits all three estimates on the fallback init if any one fails.
    # a position's three fits therefore share an init
    logliks = [ll.copy() for ll in logliks]
    if not converged.all():
        for ll, d in zip(logliks, data):
            if d.shape[1] > 1:
                ll[~converged] = fit(d[~converged], fallback=True)[1]

    return logliks[0] + logliks[1] - logliks[2], converged


def score(counts1, counts2, qual):
    """JACUSA2's call-2 score column, over ACGT counts."""
    return log_likelihood_ratio(
        pseudocounts(_with_replicates(counts1), qual),
        pseudocounts(_with_replicates(counts2), qual))


def score_indel(counts1, counts2):
    """JACUSA2's insertion_score and deletion_score columns.

    Counts are (N, 2) of reads with the indel against reads lacking it. The
    indel path adds a flat pseudo-count where the base path spreads a base
    error, and reports twice the likelihood ratio.
    """
    stat, converged = log_likelihood_ratio(
        _with_replicates(np.asarray(counts1, float) + INDEL_PSEUDOCOUNT),
        _with_replicates(np.asarray(counts2, float) + INDEL_PSEUDOCOUNT))

    return 2 * stat, converged


def scorer_for(error_type, qual):
    """The scorer for one error type. qual is ignored by the indel scorers."""
    if error_type == "mismatch":
        return partial(score, qual=qual)

    return score_indel


# ----------------------------------------------------------- the permutation
# pool the reads of every replicate at one position
# permute them back to each replicate at its own depth, without replacement
# read the observed score against the scores of the permutations
# a table, or read split, is one such assignment, (R, K) counts


def enumerate_splits(pooled, size, limit):
    """Every way to take `size` reads from the pooled counts, as (S, K).

    Returns None past limit. Categories are added one at a time. A partial
    vector is kept while it holds at most `size` reads and the categories
    still to come can make up the rest.
    """
    pooled = np.asarray(pooled, dtype=np.int64)
    fillable = np.concatenate([np.cumsum(pooled[::-1])[::-1][1:], [0]])
    rows, taken = np.zeros((1, 0), dtype=np.int64), np.zeros(1, dtype=np.int64)

    for category, cap in enumerate(pooled):
        # allocation comes before pruning, and a deep position can expand far
        # past its survivors. 8x limit is the headroom allowed
        if rows.shape[0] * (cap + 1) > 8 * limit:
            return None

        added = np.tile(np.arange(cap + 1), len(taken))
        rows = np.column_stack([np.repeat(rows, cap + 1, axis=0), added])
        taken = np.repeat(taken, cap + 1) + added
        keep = (taken <= size) & (taken + fillable[category] >= size)
        rows, taken = rows[keep], taken[keep]

    return rows if len(rows) <= limit else None


def enumerate_tables(pooled, depths, limit):
    """Every table of the pooled counts with the given group depths fixed.

    Returns (S, R, K), or None past limit tables. The last group takes what
    the others leave.
    """
    pooled = np.asarray(pooled, dtype=np.int64)
    tables = np.zeros((1, 0, len(pooled)), dtype=np.int64)
    for depth in depths[:-1]:
        grown, held = [], 0
        for table in tables:
            options = enumerate_splits(pooled - table.sum(0), depth, limit)
            if options is None or held + len(options) > limit:
                return None
            held += len(options)
            grown.append(np.concatenate(
                [np.repeat(table[None], len(options), axis=0),
                 options[:, None]], axis=1))
        tables = np.concatenate(grown)

    return np.concatenate([tables, (pooled - tables.sum(1))[:, None]], axis=1)


def table_pmf(tables, pooled):
    """Multivariate hypergeometric probability of each (S, R, K) table."""
    mass = np.ones(len(tables))
    left = np.broadcast_to(pooled, tables[:, 0].shape)
    for group in range(tables.shape[1] - 1):
        mass *= multivariate_hypergeom.pmf(x=tables[:, group], m=left,
                                           n=tables[:, group].sum(-1))
        left = left - tables[:, group]

    return mass


def exact_tail(counts, first, observed, scorer, limit):
    """Exact p for every position whose tables fit under limit.

    counts is (N, R, K), its first `first` groups condition 1. Returns the p
    and the table count. A count of zero marks a position left to the draws.
    """
    # counts arrive float for the scorer. the pmf takes integers only
    counts = counts.astype(np.int64)
    pooled, depths = counts.sum(1), counts.sum(-1)
    pvalue = np.full(len(counts), np.nan)
    n_splits = np.zeros(len(counts), dtype=np.int64)

    def flush(groups, positions):
        # reduceat reads arr[start] where two starts coincide. the observed
        # table is always in its own support, which keeps every span positive
        # and the starts strictly increasing
        spans = [len(g) for g in groups]
        starts = np.concatenate([[0], np.cumsum(spans)[:-1]])
        rows = np.concatenate(groups).astype(float)
        row_position = np.repeat(positions, spans)
        null = scorer(rows[:, :first], rows[:, first:])[0]
        mass = np.concatenate([table_pmf(g, pooled[i])
                               for g, i in zip(groups, positions)])
        tail = np.where(null >= observed[row_position] - TIE_TOLERANCE,
                        mass, 0.0)
        # the sum over a whole support is 1 give or take float noise, and a p
        # of 1 + 1e-13 is out of range for the correction downstream
        pvalue[positions] = np.clip(np.add.reduceat(tail, starts), 0.0, 1.0)

    groups, positions, held = [], [], 0
    for i in range(len(counts)):
        group = enumerate_tables(pooled[i], depths[i], limit)
        if group is None:
            continue

        groups.append(group)
        positions.append(i)
        n_splits[i] = len(group)
        held += len(group)
        if held >= BATCH_TABLES:
            flush(groups, positions)
            groups, positions, held = [], [], 0

    if groups:
        flush(groups, positions)

    return pvalue, n_splits


def permute_reads(counts, rng, draws):
    """Random permutations of (N, R, K) counts, as (N, draws, R, K) tables.

    The reads of all groups are pooled and reassigned, each group keeping its
    own depth, drawn from what the earlier groups left, one category at a time.
    """
    positions, groups, categories = counts.shape
    counts = counts.astype(np.int64)
    left = np.repeat(counts.sum(1)[:, None], draws, axis=1)
    out = np.zeros((positions, draws, groups, categories), dtype=np.int64)
    for group in range(groups - 1):
        size = np.repeat(counts[:, group].sum(-1)[:, None], draws, axis=1)
        rest = left.sum(-1)
        for category in range(categories - 1):
            rest -= left[..., category]
            taken = rng.hypergeometric(left[..., category], rest, size)
            out[:, :, group, category] = taken
            size -= taken
        out[:, :, group, -1] = size
        left -= out[:, :, group]
    out[:, :, -1] = left

    return out


def score_tables(permuted, first, scorer):
    """Score (P, D, R, K) permuted tables, fitting each distinct one once.

    Returns (P, D). A table is keyed by its position as well as its counts,
    since the same counts at another position come from a different pool.
    """
    positions, draws, groups, categories = permuted.shape
    flat = permuted.reshape(-1, groups, categories)
    key = np.column_stack([np.repeat(np.arange(positions), draws),
                           flat.reshape(len(flat), -1)])
    _, index, back = np.unique(key, axis=0, return_index=True,
                               return_inverse=True)
    tables = flat[index].astype(float)

    return scorer(tables[:, :first], tables[:, first:])[0][back].reshape(
        positions, draws)


def permutation_test(conditions, scorer, seed, draws=MIN_PERMUTATIONS):
    """Score every position against permutations of its own reads.

    conditions is the two conditions' counts, each (N, R, K), one row per
    replicate. scorer comes from scorer_for and scores observed and null.

    The p is exact where at most `draws` tables exist, and drawn elsewhere.
    n_higher counts the permutations at or above the observed score, for
    redraw to extend.
    """
    counts1, counts2 = conditions
    first = counts1.shape[1]
    counts = np.concatenate([counts1, counts2], axis=1)
    rng = np.random.default_rng(seed)
    observed, _ = scorer(counts1, counts2)
    pvalue, n_splits = exact_tail(counts, first, observed, scorer, draws)

    null = score_tables(permute_reads(counts, rng, draws), first, scorer)
    n_higher = (null >= observed[:, None]).sum(-1)
    drawn = n_splits == 0
    pvalue[drawn] = (1 + n_higher[drawn]) / (draws + 1)

    return {"score": observed,
            "score_minus_max_null": np.clip(observed - null.max(-1), 0.0,
                                            None),
            "pvalue": pvalue, "n_splits": n_splits, "p_exact": n_splits > 0,
            "n_permutations": np.where(drawn, draws, 0), "n_higher": n_higher}


def count_higher(task):
    """Permutations at or above the observed score, per position, for one batch.

    Returns the error type, the batch offset, the count per site, and the
    0-based index in the batch of each of the first `min_higher` such
    permutations, -1 where the batch holds fewer.
    """
    (counts, first, observed, error_type, qual, draws, seed, done,
     min_higher) = task
    null = score_tables(permute_reads(counts, np.random.default_rng(seed),
                                      draws), first, scorer_for(error_type, qual))
    hits = (null >= observed[:, None]).cumsum(-1)
    order = np.arange(1, min_higher + 1)
    index = (hits[:, :, None] >= order).argmax(1)

    return error_type, done, hits[:, -1], np.where(hits[:, -1:] >= order,
                                                   index, -1)


def redraw(results, conditions, args, mapper):
    """Draw more permutations for sites with too few scoring at or above them.

    A site stops at the permutation that brings the count scoring at or
    above its observed score to args.min_higher_scoring, with p = that count
    / permutations drawn. Otherwise it stops at args.max_permutations, with
    p = (count + 1) / (permutations + 1). A site with enough such
    permutations among the first draws keeps its p. Updates results in
    place. conditions maps each error type to its counts.
    """
    min_higher = args.min_higher_scoring
    for round_ in itertools.count(1):
        drawing = {t: np.flatnonzero((r["n_permutations"] > 0)
                                     & (r["n_higher"] < min_higher)
                                     & (r["n_permutations"]
                                        < args.max_permutations))
                   for t, r in results.items()}
        drawing = {t: again for t, again in drawing.items() if again.size}
        if not drawing:
            return

        started = time.monotonic()
        tasks, extra = [], {}
        for error_type, again in drawing.items():
            result = results[error_type]
            first = conditions[error_type][0].shape[1]
            counts = np.concatenate(conditions[error_type], axis=1)[again]
            total = result["n_permutations"][again[0]]
            extra[error_type] = min(total * (REDRAW - 1),
                                    args.max_permutations - total)
            # batches of at most BATCH_TABLES tables bound the memory, and a
            # fixed batch keeps the result the same at any process count
            batch = max(1, BATCH_TABLES // again.size)
            seed = [args.seed, list(ERROR_TYPES).index(error_type), round_]
            tasks += [(counts, first, result["score"][again], error_type,
                       args.qual, min(batch, extra[error_type] - done),
                       [*seed, done], done, min_higher)
                      for done in range(0, extra[error_type], batch)]
        print(f"round {round_}: {sum(map(len, drawing.values())):,} tests "
              f"with fewer than {min_higher} permutations at or above the "
              f"observed score; {len(tasks):,} "
              "batches", flush=True)

        parts = {t: [] for t in drawing}
        for error_type, *part in mapper(count_higher, tasks):
            parts[error_type].append(part)
        stopped = 0
        for error_type, again in drawing.items():
            result = results[error_type]
            total = result["n_permutations"][again[0]]
            # 1-based permutation number of each one at or above, sorted
            draw = np.sort(np.concatenate(
                [np.where(index >= 0, total + done + index + 1, np.inf)
                 for done, _, index in parts[error_type]], axis=1), axis=1)
            need = min_higher - result["n_higher"][again]
            at = draw[np.arange(again.size), need - 1]
            result["n_higher"][again] += sum(h for _, h, _ in parts[error_type])
            done_now = np.isfinite(at)
            stopped += int(done_now.sum())

            n = np.where(done_now, at, total + extra[error_type]).astype(int)
            result["n_permutations"][again] = n
            result["pvalue"][again] = np.where(
                done_now, min_higher / n,
                (1 + result["n_higher"][again]) / (n + 1))
        print(f"round {round_} done in "
              f"{(time.monotonic() - started) / 60:.1f} min, {stopped:,} tests "
              f"stopped at {min_higher} permutations at or above",
              flush=True)


# ------------------------------------------------------------ the input file


def read_header(path):
    """The JACUSA2 command line and the column names, from the two # lines.

    Raises InputError when either line is absent.
    """
    command, columns = None, None
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("##"):
                command = line[2:].strip()
            elif line.startswith("#contig"):
                columns = line[1:].rstrip("\n").split("\t")

    if columns is None:
        raise InputError(f"{path}: no '#contig' header line.")
    if command is None:
        raise InputError(f"{path}: no '##' line, where the -Q setting is.")

    return command, columns


def read_quality(command):
    """The -Q value from the command. Where -Q is absent, None."""
    quality = re.search(r"-Q\s+(\d+)", command)

    return int(quality.group(1)) if quality else None


def check_columns(columns):
    """Raise InputError unless there are two conditions and every count column.

    Conditions are read off the bases<condition><replicate> columns.
    """
    found = sorted({int(c[5]) for c in columns
                    if re.fullmatch(r"bases\d\d+", c)})
    if found != [1, 2]:
        raise InputError(f"expected two conditions, found {found or 'none'}.")

    needed = {f"{counts}{c}1" for counts in ERROR_TYPES.values()
              for c in (1, 2)}
    missing = sorted(needed - set(columns))
    if missing:
        raise InputError(f"missing columns: {', '.join(missing)}. The call-2 "
                         "run needs -D -i and -f X:reads:nonref_ratio:"
                         "insertion_ratio:deletion_ratio.")


def acgt(column):
    """An (N, 4) array of ACGT counts, from a comma-separated bases column."""
    return np.array([[int(v) for v in row.split(",")] for row in column],
                    dtype=float)


def counts_for(table, error_type):
    """Both conditions' counts for one error type, each (N, R, K).

    Mismatches read the ACGT vectors, indels the reads with the indel against
    the reads spanning the position. Each replicate is one row.
    """
    counts = ERROR_TYPES[error_type]
    conditions = []
    for condition in (1, 2):
        rows = []
        for name in table.filter(regex=rf"^bases{condition}\d+$").columns:
            bases = acgt(table[name])
            if error_type != "mismatch":
                replicate = name.removeprefix("bases")
                found = table[f"{counts}{replicate}"].to_numpy(float)
                spanning = (bases.sum(1)
                            + table[f"deletions{replicate}"].to_numpy(float))
                bases = np.c_[found, spanning - found]
            rows.append(bases)
        conditions.append(np.stack(rows, axis=1))

    return tuple(conditions)


def read_table(path, columns):
    """Every column that JACUSA2 wrote, named by the header line."""
    return pd.read_csv(path, sep="\t", comment="#", names=columns,
                       na_values=["*"])


# --------------------------------------------------------------- the command


def permute_chunk(task):
    """Score one chunk. The index comes back so results can be reordered."""
    index, first, second, draws, seed, error_type, qual = task

    return index, permutation_test((first, second),
                                   scorer_for(error_type, qual), seed,
                                   draws=draws)


def report_progress(done, total, started, error_type):
    """Print one line per tenth of the chunks. Equal chunk cost is assumed."""
    if done % max(total // 10, 1) and done != total:
        return

    elapsed = time.monotonic() - started
    left = elapsed / done * (total - done)
    print(f"{error_type}: chunk {done:,} of {total:,}, "
          f"{elapsed / 60:.1f} min in, {left / 60:.1f} min left", flush=True)


def permute(conditions, error_type, args, mapper):
    """Permutation results for one error type, chunked through mapper.

    Each chunk is seeded from its index. The process count does not change the
    numbers.
    """
    first, second = conditions
    tasks = [(n, first[i:i + CHUNK], second[i:i + CHUNK],
              args.min_permutations, args.seed + n, error_type, args.qual)
             for n, i in enumerate(range(0, len(first), CHUNK))]

    started = time.monotonic()
    results = [None] * len(tasks)
    for done, (index, result) in enumerate(mapper(permute_chunk, tasks), 1):
        results[index] = result
        report_progress(done, len(tasks), started, error_type)

    return {key: np.concatenate([chunk[key] for chunk in results])
            for key in results[0]}


def parse_args():
    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument("jacusa2_out", help="a JACUSA2 call-2 output file")
    parser.add_argument("-o", "--output", required=True, help="TSV to write")
    parser.add_argument("--seed", type=int, default=SEED,
                        help=f"seeds the permutations, default {SEED}")
    parser.add_argument("--qual", type=int,
                        help="base quality the mismatch pseudocounts "
                             "assume, taken from the -Q in the file's header "
                             "when left out")
    parser.add_argument("--min-permutations", type=int,
                        default=MIN_PERMUTATIONS,
                        help="permutations drawn for every site, default "
                             f"{MIN_PERMUTATIONS}. A site with at most this "
                             "many possible read splits gets an exact p "
                             "instead, and score_minus_max_null subtracts "
                             "the highest of this many permuted scores")
    parser.add_argument("--min-higher-scoring", type=int,
                        default=MIN_HIGHER_SCORING,
                        help="draw permutations for a site until this many "
                             "score at or above its observed score, up to "
                             "--max-permutations, default "
                             f"{MIN_HIGHER_SCORING}")
    parser.add_argument("--max-permutations", type=int,
                        default=MAX_PERMUTATIONS,
                        help=f"rounds of {REDRAW}x more permutations stop at "
                             f"this many, default {MAX_PERMUTATIONS:,}")
    parser.add_argument("--processes", type=int, default=os.cpu_count(),
                        help="default every core. the result is the same at "
                             "any number")
    return parser, parser.parse_args()


def main():
    parser, args = parse_args()
    run_started = time.monotonic()

    try:
        command, columns = read_header(args.jacusa2_out)
        check_columns(columns)
    except InputError as problem:
        sys.exit(f"{parser.prog}: {problem}")

    if args.qual is None:
        args.qual = read_quality(command)
    if args.qual is None:
        sys.exit(f"{parser.prog}: expected -Q in the header, found none. Pass "
                 "--qual with the base quality the reads carry.")

    columns[0] = "contig"
    out = read_table(args.jacusa2_out, columns)
    # one value per subsampling run is 99 per cent of the file. that statistic
    # is replaced by the permutation FDR. without a fixed seed, its columns
    # differ between runs over the same data
    dropped = [c for c in out.columns if c.endswith("_subsampled")]
    out = out.drop(columns=dropped)

    print(f"{len(out):,} positions, {args.min_permutations} to "
          f"{args.max_permutations:,} permutations, -Q {args.qual}, "
          f"{args.processes} process(es)", flush=True)
    if dropped:
        print(f"dropped {', '.join(dropped)}", flush=True)

    conditions = {t: counts_for(out, t) for t in ERROR_TYPES}
    with Pool(args.processes) if args.processes > 1 else nullcontext() as pool:
        mapper = pool.imap_unordered if pool else map
        results = {t: permute(conditions[t], t, args, mapper)
                   for t in ERROR_TYPES}
        redraw(results, conditions, args, mapper)
    # one BH family over every error type. a position is called on its
    # lowest FDR from this one correction
    pvalues = np.concatenate([r["pvalue"] for r in results.values()])
    rates = np.split(false_discovery_control(pvalues, method="bh"),
                     len(results))
    for (error_type, result), fdr in zip(results.items(), rates):
        result["fdr"] = fdr

        added = {f"perm_{error_type}_{name}": result[name]
                 for name in ADDED_COLUMNS}
        clash = sorted(set(added) & set(out.columns))
        if clash:
            sys.exit(f"{parser.prog}: {', '.join(clash)} already in the input.")

        out = out.assign(**added)
        exact = int(result["p_exact"].sum())
        print(f"{error_type}: {exact:,} of {len(out):,} p enumerated",
              flush=True)

    out.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} in {(time.monotonic() - run_started) / 60:.1f} "
          "min", flush=True)


if __name__ == "__main__":
    main()
