"""plot_norm_heatmap.py - a heatmap from the permutation output alone

input:   the tsv that add_permutation_scores.py writes, one row per position
output:  --out, a pdf with tRNAs down the page and positions across, and
         the same page as a png beside it for the report
usage:   python plot_norm_heatmap.py WT_KO.perm.tsv --out heatmap.pdf
env:     python3 with numpy, pandas, matplotlib
notes:   a position is called on its lowest FDR. a tile's colour is the
         observed score minus the highest permutation score, summed over
         the three error types. its size says whether the FDR is below --fdr. a position that the tRNA lacks is drawn grey. --sprinzl
         renames the positions, --trnas keeps a subset, --mark draws a red box
"""
import os
import re
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.textpath import TextPath
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, MaxNLocator

TYPES = ("mismatch", "deletion", "insertion")
# ggplot's yellow to blue in Lab, sampled at nine steps
CMAP = LinearSegmentedColormap.from_list(
    "yellow_blue", ["#FFFF00", "#F5E153", "#EAC476", "#DBA792", "#CA8AAA",
                    "#B46DC0", "#9850D5", "#7030EA", "#0000FF"])
EDGE, MISSING, AXIS, GRID = "#D8C180", "#BEBEBE", "#4D4D4D", "#EBEBEB"
# a cell the tRNA has no position for carries no score, and comes out grey
CMAP.set_bad(MISSING)
# every length below is a point, measured off a plot_heatmap.R page
PITCH = 7.5
CALLED, UNCALLED = 7.13 / PITCH, 2.25 / PITCH
PAD_LEFT, PAD_BOTTOM = 17, 18.9
PAD_TOP, PAD_RIGHT, PAD_TITLE, SIDE = 62.5, 0, 17, 10
# the coverage panel: its width, the gap that holds the second tRNA labels,
# one condition's bar, and the two colours plot_heatmap.R fills them with
COVER_W, COVER_GAP, BAR_PT = 101, 6, 2.94
CONDITION = ("#BF568B", "#BCD357")
TICK, AXIS_TITLE, TITLE = 7, 9, 11
# the key band, as offsets from the left of the guide group
BAND_W, BAR_X, BAR_W, BAR_Y, BAR_H = 264, 27, 90, 25, 12
SMALL_X, KEY_GAP = 144, 13
KEY_Y, KEY_SMALL, KEY_LARGE = 28.6, 3.35, 8.08


def coverage(table, condition):
    """Base calls plus deletions in one condition, summed over replicates."""
    bases = table.filter(regex=rf"^bases{condition}\d+$").stack()
    bases = bases.str.split(",", expand=True).astype(float).sum(axis=1)
    deletions = table.filter(regex=rf"^deletions{condition}\d+$").sum(axis=1)
    return bases.groupby(level=0).sum() + deletions


def scored(table, fdr):
    """Add the normalized score and whether the position is called."""
    rates = table[[f"perm_{t}_fdr" for t in TYPES]]
    margins = table[[f"perm_{t}_score_minus_max_null" for t in TYPES]]
    return table.assign(norm_score=margins.sum(axis=1),
                        called=rates.min(axis=1) < fdr,
                        position=table.start + 1,
                        coverage1=coverage(table, 1),
                        coverage2=coverage(table, 2))


def relabel(table, sprinzl):
    """Carry the Sprinzl label over, keeping the positions that have one."""
    if not sprinzl:
        return table, "Position"

    mapping = pd.read_csv(sprinzl, sep="\t", usecols=["trna", "seq_position",
                                                      "sprinzl"])
    mapping = mapping.rename(columns={"trna": "contig",
                                      "seq_position": "position"})
    table = table.merge(mapping, on=["contig", "position"], how="left")
    table = table.dropna(subset=["sprinzl"])
    return table.assign(position=table.sprinzl), "Sprinzl pos."


def keep(table, args):
    """Apply the tRNA, depth and label filters, and name the x axis."""
    table, xlabel = relabel(table, args.sprinzl)
    if args.trnas:
        table = table[table.contig.str.contains(args.trnas, regex=True)]
    # the same rule plot_heatmap.R applies: one position clears it in both
    deep = table[(table.reads11 > args.min_reads)
                 & (table.reads21 > args.min_reads)].contig.unique()
    table = table[table.contig.isin(deep)]
    if args.remove_prefix:
        table = table.assign(
            contig=table.contig.str.removeprefix(args.remove_prefix))
    # an unlabelled position is a dot, and takes no column of its own
    table = table.assign(position=table.position.astype(str))
    return table[table.position != "."], xlabel


def ordered(table):
    """Position labels in their order along the tRNA."""
    order = table.groupby("position", sort=False).start.median()
    return [str(name) for name in order.sort_values().index]


def grid(table, trnas, positions):
    """One row per cell, with its place on the axes. A missing cell is NaN."""
    every = pd.MultiIndex.from_product([trnas, positions],
                                       names=["contig", "position"])
    cells = (table.set_index(["contig", "position"])
                  .reindex(every).reset_index())
    return cells.assign(x=cells.position.map(positions.index) + .5,
                        y=cells.contig.map(trnas.index) + .5)


def text_width(strings, size):
    """The widest of these strings, in points."""
    return max(TextPath((0, 0), string, size=size).get_extents().width
               for string in strings)


def panel(columns, rows, names, labels, *, coverage):
    """A figure whose cells are PITCH points square, with R's margins.

    Returns the figure, the heatmap axes, and the coverage axes or None.
    """
    name_w = text_width(names, TICK)
    left = name_w + PAD_LEFT + SIDE
    bottom = text_width(labels, TICK) + PAD_BOTTOM
    aside = name_w + COVER_GAP + COVER_W if coverage else 0
    width = left + columns * PITCH + aside + PAD_RIGHT + SIDE
    height = PAD_TOP + rows * PITCH + bottom

    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes([left / width, bottom / height,
                       columns * PITCH / width, rows * PITCH / height])
    if not coverage:
        return fig, ax, None

    bars = fig.add_axes([(width - SIDE - PAD_RIGHT - COVER_W) / width,
                         bottom / height, COVER_W / width,
                         rows * PITCH / height])
    return fig, ax, bars


def thousands(value, _):
    return f"{value / 1000:g}K" if value >= 1000 else f"{value:g}"


def coverage_bars(fig, bars, table, trnas, conditions):
    """Median coverage per tRNA, one bar per condition, beside the heatmap."""
    medians = table.groupby("contig")[["coverage1", "coverage2"]].median()
    medians = medians.reindex(trnas)
    height = BAR_PT / PITCH
    for offset, column, colour in ((-height / 2, "coverage1", CONDITION[0]),
                                   (height / 2, "coverage2", CONDITION[1])):
        bars.barh([index + .5 + offset for index in range(len(trnas))],
                  medians[column], height=height, color=colour, zorder=2)

    # ggplot's grid: a major line per row and per x break, minors between
    bars.set_axisbelow(True)
    bars.grid(which="major", color=GRID, linewidth=.87)
    bars.grid(axis="x", which="minor", color=GRID, linewidth=.44)
    bars.xaxis.set_minor_locator(AutoMinorLocator(2))
    bars.set_yticks(range(len(trnas) + 1), minor=True)
    # the minor breaks carry a grid line and no tick of their own
    bars.tick_params(which="minor", length=0)
    bars.set_ylim(len(trnas), 0)
    bars.set_yticks([index + .5 for index in range(len(trnas))], trnas,
                    fontsize=TICK)
    # the names repeat in the gap, between the heatmap and the bars
    bars.tick_params(axis="y", length=2.2, width=.87, color="#333333",
                     labelsize=TICK, pad=2.5)
    bars.set_xlabel("median coverage", fontsize=AXIS_TITLE)
    bars.xaxis.set_major_formatter(FuncFormatter(thousands))
    bars.tick_params(axis="x", labelsize=TICK, length=2.2, width=.87,
                     color="#333333", labelcolor=AXIS, pad=1.5)
    for side in ("top", "right", "left", "bottom"):
        bars.spines[side].set_color("#333333")
        bars.spines[side].set_linewidth(.5)

    keys = [Rectangle((0, 0), 1, 1, facecolor=colour, edgecolor="none")
            for colour in CONDITION]
    bars.legend(keys, conditions, frameon=False, fontsize=TICK,
                loc="lower left", bbox_to_anchor=(.18, 1.02),
                handlelength=1.2, handleheight=1.2, labelspacing=.45,
                handletextpad=.5)


def key_band(fig, ax, scale, fdr):
    """The score bar and the two tile sizes, centred over the panel."""
    width, height = fig.get_figwidth() * 72, fig.get_figheight() * 72
    box = ax.get_position()
    group = (box.x0 + box.x1) / 2 * width - BAND_W / 2

    def spot(x, y, w, h):
        return [x / width, 1 - (y + h) / height, w / width, h / height]

    cax = fig.add_axes(spot(group + BAR_X, BAR_Y, BAR_W, BAR_H))
    bar = fig.colorbar(ScalarMappable(scale, CMAP), cax=cax,
                       orientation="horizontal",
                       ticks=MaxNLocator(nbins=6, steps=[1, 2, 5, 10]))
    # the breaks cut into the bar itself, as ggplot's colourbar draws them
    bar.ax.tick_params(labelsize=TICK, direction="in", length=3, width=.8,
                       color="white", labelcolor="black", pad=2,
                       top=True, bottom=True)
    bar.outline.set_visible(False)

    def text(x, y, string, **kwargs):
        fig.text(x / width, 1 - y / height, string, fontsize=TICK, **kwargs)

    text(group + BAR_X - 3, BAR_Y + BAR_H / 2, "score", ha="right",
         va="center")
    # each key is a square, its label, then a gap the label's own width sets
    at = group + SMALL_X
    for side, sign in ((KEY_SMALL, ">="), (KEY_LARGE, "<")):
        fig.add_artist(Rectangle(
            spot(at - side / 2, KEY_Y - side / 2, side, side)[:2],
            side / width, side / height, facecolor="none", edgecolor=AXIS,
            linewidth=.7))
        label = f"FDR {sign} {fdr}"
        text(at + KEY_GAP, KEY_Y + 2.4, label)
        at += KEY_GAP + text_width([label], TICK) + KEY_GAP


def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("perm")
    parser.add_argument("--out", required=True)
    parser.add_argument("--fdr", type=float, default=.05)
    parser.add_argument("--min_reads", type=float, default=10,
                        help="a tRNA needs one position above this in both "
                             "conditions")
    parser.add_argument("--trnas", help="regex of tRNAs to keep")
    parser.add_argument("--mark", help="regex over <trna>:<position>")
    parser.add_argument("--sprinzl",
                        help="a scores_sprinzl.tsv of the same contrast")
    parser.add_argument("--remove_prefix", default="",
                        help="dropped from every tRNA name, e.g. tRNA-")
    parser.add_argument("--title", default="")
    parser.add_argument("--conditions",
                        help="two names, comma separated, which add the "
                             "median coverage panel, e.g. WT,QTRT1")
    args = parser.parse_args()
    if args.conditions:
        args.conditions = args.conditions.split(",")

    table, xlabel = keep(scored(pd.read_csv(args.perm, sep="\t"), args.fdr),
                         args)
    if table.empty:
        parser.error("no tRNA survives the filters")

    trnas = sorted(table.contig.unique())
    positions = ordered(table)
    cells = grid(table, trnas, positions)

    scale = Normalize(0, table.norm_score.max() or 1)
    fig, ax, bars = panel(len(positions), len(trnas), trnas, positions,
                          coverage=bool(args.conditions))
    called = cells.called.fillna(value=False)
    side = np.where(called, CALLED, UNCALLED) * PITCH
    ax.scatter(cells.x, cells.y, s=side ** 2, marker="s",
               c=CMAP(scale(np.ma.masked_invalid(cells.norm_score))),
               edgecolors=EDGE, linewidths=.25)
    if args.mark:
        marked = cells[cells.apply(
            lambda cell: bool(re.search(args.mark,
                                        f"{cell.contig}:{cell.position}")),
            axis=1)]
        # the box takes the tile's own size, as plot_heatmap.R draws it
        ax.scatter(marked.x, marked.y, marker="s", facecolors="none",
                   s=side[marked.index] ** 2, edgecolors="red", linewidths=.5)
    ax.set_xlim(0, len(positions)), ax.set_ylim(len(trnas), 0)
    ax.set_xticks([index + .5 for index in range(len(positions))], positions,
                  rotation=90, fontsize=TICK)
    ax.set_yticks([index + .5 for index in range(len(trnas))], trnas,
                  fontsize=TICK)
    ax.set_xlabel(xlabel, fontsize=AXIS_TITLE)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=1.6, width=.5, colors=AXIS, pad=1.5)
    if args.title:
        fig.text((ax.get_position().x0 * fig.get_figwidth() * 72 - PAD_TITLE)
                 / (fig.get_figwidth() * 72),
                 1 - 11 / (fig.get_figheight() * 72), args.title,
                 fontsize=TITLE, va="top")

    if bars is not None:
        coverage_bars(fig, bars, table, trnas, args.conditions)
    key_band(fig, ax, scale, args.fdr)
    # the page is already the size the margins define, and tight would resize it
    fig.savefig(args.out)
    fig.savefig(os.path.splitext(args.out)[0] + ".png", dpi=150)
    called = int(table.called.sum())
    print(f"wrote {args.out}: {len(trnas)} tRNAs, {len(positions)} positions, "
          f"{called:,} of {len(table):,} called at FDR < {args.fdr}")


if __name__ == "__main__":
    main()
