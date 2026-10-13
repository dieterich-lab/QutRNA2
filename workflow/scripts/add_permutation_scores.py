"""add_permutation_scores.py - permutation p and BH q over a JACUSA2 call-2 run

input:   a JACUSA2 call-2 output file, two conditions of one replicate each
output:  a TSV of the JACUSA2 columns bar the subsampled ones, plus the
         permutation columns
usage:   python add_permutation_scores.py JACUSA2.out -o perm.tsv --seed 42
env:     python3 with numpy, scipy, pandas. JACUSA2 is not called
notes:   new columns are perm_<type>_<name>. type is mismatch, deletion or
         insertion. name is score, score_minus_max_null, pvalue, qvalue,
         p_exact or n_splits
"""
import os
import re
import sys
import time
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.special import digamma, gammaln, polygamma
from scipy.stats import (
    false_discovery_control,
    multivariate_hypergeom,
    pearsonr,
)

CATEGORIES = 4
ESTIMATED_ERROR = 0.01
MAX_Q = 40
MAX_ITERATIONS = 100
EPSILON = 0.001
BACKTRACK_LAMBDA = 0.9
INDEL_PSEUDOCOUNT = 1.0

# past this many re-splits a position is sampled
SPLIT_LIMIT = 1_000_000
# one fit call over many positions pays the setup once
BATCH_TABLES = 250_000
# the slack keeps the observed split inside its own tail
TIE_TOLERANCE = 1e-9
POSITION_CHUNK = 250
# 25 keeps a short file above one chunk per process
CHUNK = 25
DRAWS = 1000
MIN_R = 0.99

# JACUSA2's own score column, and the prefix of its count columns
ERROR_TYPES = {
    "mismatch": ("score", "bases"),
    "deletion": ("deletion_score", "deletions"),
    "insertion": ("insertion_score", "insertions"),
}
ADDED_COLUMNS = ["score", "score_minus_max_null", "pvalue", "qvalue",
                 "p_exact", "n_splits"]


class InputError(Exception):
    """The input file cannot be scored, with the reason in the message."""


# ---------------------------------------------------------------- the score
# JACUSA2's Dirichlet-multinomial call-2 score, following
# lib/estimate/MinkaEstimateDirMultAlpha.java over all positions at once


# the cap matches the Java the port is diffed against
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

    Returns the proposed alpha, and a mask marking the positions that lack an
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

    logliks, fit_converged = zip(*(fit(d)[1:] for d in data))
    converged = np.logical_and.reduce(fit_converged)

    # JACUSA2 refits all three estimates on the fallback init if any one fails.
    # a position's three fits therefore share an init
    logliks = [ll.copy() for ll in logliks]
    if not converged.all():
        for ll, d in zip(logliks, data):
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
# pool the counts of one position
# re-split into groups of the original sizes, without replacement
# read the observed score against the re-split scores


def enumerate_splits(pooled, size, limit=SPLIT_LIMIT):
    """Every distinct re-split of one position, as (S, K). Past limit, None.

    Categories are added one at a time. A partial vector survives while it has
    spent at most the group and the rest can still fill it, so every surviving
    prefix completes.
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


def exact_tail(pooled, counts1, observed, scorer, limit=SPLIT_LIMIT):
    """Exact p for every position whose re-splits fit under limit.

    Returns the p and the re-split count. A count of zero marks a position left
    to the sampler.
    """
    pvalue = np.full(len(pooled), np.nan)
    n_splits = np.zeros(len(pooled), dtype=np.int64)

    def flush(groups, positions):
        # reduceat reads arr[start] where two starts coincide. the observed
        # split is always in its own support, which keeps every span positive
        # and the starts strictly increasing
        spans = [len(g) for g in groups]
        starts = np.concatenate([[0], np.cumsum(spans)[:-1]])
        rows = np.concatenate(groups)
        row_position = np.repeat(positions, spans)
        null = scorer(rows.astype(float),
                      (pooled[row_position] - rows).astype(float))[0]
        mass = multivariate_hypergeom.pmf(x=rows, m=pooled[row_position],
                                          n=rows.sum(-1))
        tail = np.where(null >= observed[row_position] - TIE_TOLERANCE,
                        mass, 0.0)
        # the sum over a whole support is 1 give or take float noise, and a p
        # of 1 + 1e-13 is out of range for the correction downstream
        pvalue[positions] = np.clip(np.add.reduceat(tail, starts), 0.0, 1.0)

    groups, positions, held = [], [], 0
    for i in range(len(pooled)):
        group = enumerate_splits(pooled[i], int(counts1[i].sum()), limit)
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


def split_once(pooled, sizes, rng, draws):
    """Re-split pooled counts into two groups of the original sizes.

    pooled is (N, 4) and sizes is (N,), the read count of condition 1. Returns
    two (N, draws, 4) arrays that add back to pooled.
    """
    first = np.stack([rng.multivariate_hypergeometric(p, n, size=draws)
                      for p, n in zip(pooled.astype(int), sizes)])

    return first.astype(float), pooled[:, None, :] - first


def score_splits(left, right, scorer):
    """Score every re-split, fitting each distinct table once. (positions, draws).

    A table is keyed by its position as well as its counts, since the same
    counts at another position come from a different pool.
    """
    positions, draws, categories = left.shape
    flat_left = left.reshape(-1, categories)
    flat_right = right.reshape(-1, categories)
    key = np.column_stack([np.repeat(np.arange(positions), draws),
                           flat_left.astype(np.int64)])
    _, first, back = np.unique(key, axis=0, return_index=True,
                               return_inverse=True)
    scores = scorer(flat_left[first], flat_right[first])[0]

    return scores[back].reshape(positions, draws)


def permutation_test(conditions, qual, seed, draws=DRAWS,
                     error_type="mismatch"):
    """Score every position against the permutation null built from it.

    conditions is the two conditions' counts, each (N, K). The mismatch
    pseudocounts are set by qual, over observed and null alike.

    The p is enumerated where the re-splits fit under SPLIT_LIMIT and drawn
    where they do not, bounded at 1 / (draws + 1). p_exact marks which. The
    score always comes from the draws, so every null holds the same count.
    """
    counts1, counts2 = conditions
    scorer = scorer_for(error_type, qual)
    rng = np.random.default_rng(seed)
    observed, _ = scorer(counts1, counts2)

    pooled = counts1 + counts2
    sizes = counts1.sum(-1).astype(int)
    pvalue, n_splits = exact_tail(pooled, counts1, observed, scorer)

    null_max = np.full(len(pooled), np.nan)
    for start in range(0, len(pooled), POSITION_CHUNK):
        block = slice(start, start + POSITION_CHUNK)
        left, right = split_once(pooled[block], sizes[block], rng, draws)
        null = score_splits(left, right, scorer)
        null_max[block] = null.max(-1)

        over = n_splits[block] == 0
        if over.any():
            drawn = (1 + (null >= observed[block][:, None]).sum(-1)) / (draws + 1)
            pvalue[block] = np.where(over, drawn, pvalue[block])

    return {"score": observed,
            "score_minus_max_null": np.clip(observed - null_max, 0.0, None),
            "pvalue": pvalue, "n_splits": n_splits, "p_exact": n_splits > 0}


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


def count_replicates(columns):
    """Replicates per condition, as {condition: count}, read off bases<c><r>."""
    replicates = {}
    for column in columns:
        found = re.fullmatch(r"bases(\d)(\d+)", column)
        if found:
            condition = int(found.group(1))
            replicates[condition] = replicates.get(condition, 0) + 1

    return replicates


def check_shape(columns):
    """Raise InputError unless the result is two conditions of one replicate."""
    replicates = count_replicates(columns)
    if sorted(replicates) != [1, 2]:
        raise InputError(f"expected two conditions, found "
                         f"{sorted(replicates) or 'none'}.")

    crowded = {c: n for c, n in replicates.items() if n > 1}
    if crowded:
        listed = ", ".join(f"condition {c} has {n}" for c, n in crowded.items())
        raise InputError(f"expected one replicate per condition, {listed}. Run "
                         "JACUSA2 per replicate pair, or merge the BAMs "
                         "first.")


def check_columns(columns, error_types):
    """Raise InputError naming the -f option when a needed column is absent."""
    needed = {"bases11", "bases21"}
    for error_type in error_types:
        column, counts = ERROR_TYPES[error_type]
        needed |= {column, f"{counts}11", f"{counts}21"}
        if error_type != "mismatch":
            needed |= {"deletions11", "deletions21"}

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
    """Both conditions' counts for one error type, each (N, K).

    Mismatches read the ACGT vectors, indels the reads with the indel against
    the reads spanning the position.
    """
    _, counts = ERROR_TYPES[error_type]
    if error_type == "mismatch":
        return tuple(acgt(table[f"bases{c}1"]) for c in (1, 2))

    spanning = [acgt(table[f"bases{c}1"]).sum(1) + table[f"deletions{c}1"].values
                for c in (1, 2)]

    return tuple(np.c_[table[f"{counts}{c}1"].values,
                       total - table[f"{counts}{c}1"].values].astype(float)
                 for c, total in zip((1, 2), spanning))


def read_table(path, columns):
    """Every column that JACUSA2 wrote, named by the header line."""
    return pd.read_csv(path, sep="\t", comment="#", names=columns,
                       na_values=["*"])


def agreement(recomputed, reported):
    """Pearson r between the recomputed score and JACUSA2's own column."""
    usable = np.isfinite(recomputed) & np.isfinite(reported)
    if usable.sum() < 2:
        return float("nan")

    return pearsonr(reported[usable], recomputed[usable])[0]


# --------------------------------------------------------------- the command


def permute_chunk(task):
    """Score one chunk. The index comes back so results can be reordered."""
    index, first, second, draws, seed, error_type, qual = task

    return index, permutation_test((first, second), qual, seed, draws=draws,
                                   error_type=error_type)


def report_progress(done, total, started, error_type):
    """Print one line per tenth of the chunks. Equal chunk cost is assumed."""
    if done % max(total // 10, 1) and done != total:
        return

    elapsed = time.monotonic() - started
    left = elapsed / done * (total - done)
    print(f"{error_type}: chunk {done:,} of {total:,}, "
          f"{elapsed / 60:.1f} min in, {left / 60:.1f} min left", flush=True)


def permute(conditions, error_type, args):
    """Permutation results for one error type, chunked across a pool.

    Each chunk is seeded from its index. The process count does not change the
    numbers.
    """
    first, second = conditions
    tasks = [(n, first[i:i + args.chunk], second[i:i + args.chunk], args.draws,
              args.seed + n, error_type, args.qual)
             for n, i in enumerate(range(0, len(first), args.chunk))]

    started = time.monotonic()
    done = 0
    if args.processes == 1:
        results = []
        for task in tasks:
            results.append(permute_chunk(task)[1])
            done += 1
            report_progress(done, len(tasks), started, error_type)
    else:
        results = [None] * len(tasks)
        with Pool(args.processes) as pool:
            for index, result in pool.imap_unordered(permute_chunk, tasks):
                results[index] = result
                done += 1
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
    parser.add_argument("--seed", type=int, required=True,
                        help="seeds the drawn re-splits")
    parser.add_argument("--qual", type=int,
                        help="base quality the mismatch pseudocounts "
                             "assume, taken from the -Q in the file's header "
                             "when left out")
    parser.add_argument("--error-types", default="mismatch,deletion,insertion",
                        help="comma separated, default all three")
    parser.add_argument("--draws", type=int, default=DRAWS,
                        help=f"re-splits drawn per position, default {DRAWS}")
    parser.add_argument("--chunk", type=int, default=CHUNK,
                        help=f"positions per task, default {CHUNK}")
    parser.add_argument("--processes", type=int, default=os.cpu_count(),
                        help="default every core. the result is the same at "
                             "any number")
    parser.add_argument("--min-r", type=float, default=MIN_R,
                        help=f"stop when the mismatch r falls below this, "
                             f"default {MIN_R}")
    args = parser.parse_args()

    args.error_types = [t.strip() for t in args.error_types.split(",")
                        if t.strip()]
    unknown = sorted(set(args.error_types) - set(ERROR_TYPES))
    if unknown:
        parser.error(f"unknown error type(s): {', '.join(unknown)}. "
                     f"Choose from {', '.join(ERROR_TYPES)}.")

    return parser, args


def main():
    parser, args = parse_args()
    run_started = time.monotonic()

    try:
        command, columns = read_header(args.jacusa2_out)
        check_shape(columns)
        check_columns(columns, args.error_types)
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
    # is replaced by the permutation q. without a fixed seed, its columns
    # differ between runs over the same data
    dropped = [c for c in out.columns if c.endswith("_subsampled")]
    out = out.drop(columns=dropped)

    print(f"{len(out):,} positions, {args.draws} draws, -Q {args.qual}, "
          f"{len(args.error_types)} error type(s), "
          f"{args.processes} process(es)",
          flush=True)
    if dropped:
        print(f"dropped {', '.join(dropped)}", flush=True)

    for error_type in args.error_types:
        column, _ = ERROR_TYPES[error_type]
        result = permute(counts_for(out, error_type), error_type, args)
        result["qvalue"] = false_discovery_control(result["pvalue"],
                                                   method="bh")

        added = {f"perm_{error_type}_{name}": result[name]
                 for name in ADDED_COLUMNS}
        clash = sorted(set(added) & set(out.columns))
        if clash:
            sys.exit(f"{parser.prog}: {', '.join(clash)} already in the input.")

        out = out.assign(**added)
        r = agreement(result["score"], out[column].to_numpy(float))
        exact = int(result["p_exact"].sum())
        print(f"{error_type}: r {r:.4f} against {column}, "
              f"{exact:,} of {len(out):,} p enumerated", flush=True)
        if error_type == "mismatch" and r < args.min_r:
            sys.exit(f"{parser.prog}: mismatch r {r:.4f} is below "
                     f"{args.min_r}. The recomputed score departs from "
                     "JACUSA2's by more than this run allows.")

    out.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {args.output} in {(time.monotonic() - run_started) / 60:.1f} "
          "min", flush=True)


if __name__ == "__main__":
    main()
