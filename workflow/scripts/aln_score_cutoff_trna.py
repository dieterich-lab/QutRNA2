#!/usr/bin/env python3
"""
compute per-trna alignment score cutoffs using one global precision target.

inputs are binned count tables from:
  bam_utils.py count-trna-wise-tag -t AS -c alignment_score

output is a tsv consumed by bam_utils.py filter --cutoffs with columns:
  trna    cutoff

notes:
- `trna` is an exact reference name from the score tables.
- this script does not build regex groups; every observed trna gets its own cutoff.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import precision_recall_curve


def load_tsv_expand(tsv_file: Path, category: str) -> pd.DataFrame:
    """Load count table and expand to one row per alignment for PR computation."""
    try:
        # Input files already contain a header: trna, alignment_score, count
        df = pd.read_csv(tsv_file, sep="\t")
        required = {"trna", "alignment_score", "count"}
        if not required.issubset(df.columns):
            missing = sorted(required - set(df.columns))
            raise ValueError(f"Missing required column(s): {', '.join(missing)}")

        # Expand weighted counts into per-alignment rows to keep PR computation explicit.
        expanded = df.loc[df.index.repeat(df["count"]), ["trna", "alignment_score"]].copy()
        expanded["category"] = category
        return expanded
    except Exception as exc:
        sys.stderr.write(f"Error loading {tsv_file}: {exc}\n")
        sys.exit(1)


def compute_cutoff(group_scores: pd.DataFrame, target_precision: float) -> float:
    """Return the lowest alignment score reaching target precision."""
    # Precision-recall is computed on binary labels where POS are real alignments.
    y_true = (group_scores["category"] == "POS").astype(int).to_numpy()
    y_score = group_scores["alignment_score"].to_numpy()

    prec, _, thresh = precision_recall_curve(y_true, y_score)

    # precision_recall_curve returns len(prec) == len(thresh) + 1.
    # Compare against prec[:-1] so indices align to threshold values.
    mask = prec[:-1] >= target_precision
    if not np.any(mask):
        # If target precision is unreachable, fall back to strictest observed score.
        cutoff = float(group_scores["alignment_score"].max())
        sys.stderr.write(
            f"Warning: cannot reach precision {target_precision:.3f}; using max score {cutoff:.1f}\n"
        )
        return cutoff

    idx = int(np.argmax(mask))
    return float(thresh[idx])




def plot_groups(plot_path: Path, plot_rows: list[pd.DataFrame]) -> None:
    """render one panel per exact trna showing real/random score densities and cutoff."""
    if not plot_rows:
        sys.stderr.write("no data to plot\n")
        return

    plot_df = pd.concat(plot_rows, ignore_index=True)
    trna_names = list(plot_df["group"].unique())
    n_groups = len(trna_names)

    # use one subplot per trna so each threshold line is easy to inspect.
    _, axes = plt.subplots(n_groups, 1, figsize=(10, 3 * n_groups), sharex=True, sharey=True)
    if n_groups == 1:
        axes = [axes]

    for ax, (trna_name, group_df) in zip(axes, plot_df.groupby("group", sort=False)):
        cutoff = float(group_df["cutoff"].iloc[0])

        # overlay real vs random score distributions for this exact trna.
        for cat, color in (("POS", "orange"), ("NEG", "blue")):
            vals = group_df.loc[group_df["category"] == cat, "alignment_score"].dropna()
            if len(vals) > 0:
                vals_df = pd.DataFrame({"alignment_score": vals})
                sns.histplot(
                    data=vals_df,
                    x="alignment_score",
                    bins=20,
                    alpha=0.65,
                    stat="density",
                    ax=ax,
                    color=color,
                    label=cat,
                )

        ax.axvline(cutoff, color="red", linestyle="--", linewidth=2, label=f"cutoff {cutoff:.1f}")
        ax.set_title(trna_name)
        ax.set_ylabel("density")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right")

    plt.xlabel("alignment score")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    sys.stderr.write(f"score plot -> {plot_path}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="per-trna alignment score cutoffs from one global precision target."
    )
    parser.add_argument("--real", required=True, type=Path, help="pos table: trna/alignment_score/count")
    parser.add_argument("--random", required=True, type=Path, help="neg table: trna/alignment_score/count")
    parser.add_argument("--precision", type=float, default=0.999, help="global precision target (0-1)")
    parser.add_argument("--min-samples", type=int, default=20, help="skip trnas with fewer samples")
    parser.add_argument("--output", required=True, type=Path, help="output cutoff tsv (trna, cutoff)")
    parser.add_argument("--score-plot", type=Path, default=None, help="optional pdf plot")
    args = parser.parse_args()

    if not (0.0 < args.precision <= 1.0):
        sys.stderr.write("error: --precision must be in (0, 1].\n")
        sys.exit(1)

    # combine real and random scores with explicit labels used by pr computation.
    real_df = load_tsv_expand(args.real, "POS")
    rand_df = load_tsv_expand(args.random, "NEG")

    # random alignments are reversed sequences, so reference names carry a '_rev' suffix.
    # strip it so NEG and POS rows from the same trna fall into the same group.
    rand_df["trna"] = rand_df["trna"].str.replace(r"_rev$", "", regex=True)

    scores = pd.concat([real_df, rand_df], ignore_index=True)

    # use all unique observed trna names directly, with one cutoff per exact name.
    unique_trnas = sorted(scores["trna"].dropna().unique().tolist())

    cutoff_rows: list[dict[str, float | str]] = []
    plot_rows: list[pd.DataFrame] = []

    for trna_name in unique_trnas:
        group_scores = scores[scores["trna"] == trna_name].copy()
        total_samples = len(group_scores)

        if total_samples < args.min_samples:
            sys.stderr.write(f"skipping {trna_name}: {total_samples} samples (< {args.min_samples})\n")
            continue

        cutoff = compute_cutoff(group_scores, args.precision)
        cutoff_rows.append({"trna": trna_name, "cutoff": cutoff})
        sys.stderr.write(f"{trna_name}: cutoff={cutoff:.1f} (N={total_samples})\n")

        if args.score_plot is not None:
            # keep compact plotting fields so the plot phase stays simple.
            group_scores["group"] = trna_name
            group_scores["cutoff"] = cutoff
            plot_rows.append(group_scores[["alignment_score", "category", "group", "cutoff"]])

    if not cutoff_rows:
        sys.stderr.write("error: no cutoffs computed\n")
        sys.exit(1)

    # keep output schema compatible with bam_utils.py filter loader.
    pd.DataFrame(cutoff_rows).to_csv(args.output, sep="\t", index=False)
    sys.stderr.write(f"cutoffs -> {args.output}\n")

    if args.score_plot is not None:
        plot_groups(args.score_plot, plot_rows)


if __name__ == "__main__":
    main()
