import csv
import os
import logging
import re
from collections import defaultdict

import click
import pandas as pd
import pysam
from Bio import AlignIO


logger = logging.getLogger(__name__)


@click.group()
def cli():
    pass


class FastaRecord:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq


def write_fasta_records(records, fname):
    with open(fname, "w") as f:
        for record in records.values():
            f.write(f">{record.id}\n")
            f.write(f"{record.seq}\n")


def _label_sort_key(label):
    """sort key placing every label in true 5'->3' Sprinzl order.
    plain numeric labels (with optional letter-suffix overflow, e.g. '60A')
    sort by (number, suffix length, suffix): suffix length before suffix
    itself so a single-letter overflow (Z) sorts before the two-letter
    overflow that follows it (AA), which plain string comparison gets wrong.
    variable-arm 'e' labels (class-ii: Leu, Ser) sit strictly between 45 and
    46: e11-e17 (5' stem) first, then e1-e5 (loop), then e21-e27 (3' stem),
    but the 3' stem is numbered in REVERSE as you read it 5'->3' (e27
    comes before e21), so its sort order needs the digit negated or e27
    would wrongly sort after e21."""
    m = re.match(r"e(\d)(\d)?([A-Za-z]*)$", label)
    if m:
        d1, d2, suffix = m.group(1), m.group(2), m.group(3)
        if d2 is None:
            rank, order = 1, int(d1)      # e1..e5: loop
        elif d1 == "1":
            rank, order = 0, int(d2)      # e11..e17: 5' stem
        else:
            rank, order = 2, -int(d2)     # e21..e27: 3' stem, reverse order
        return (45, 1, rank, order, len(suffix), suffix)
    m = re.match(r"(\d+)([A-Za-z]*)", label)
    num, suffix = m.group(1), m.group(2)
    return (int(num), 0, 0, 0, len(suffix), suffix)


def _column_has_residue(align, col_idx):
    residue_count = 0
    gap_count = 0
    other_symbols = []
    unusual_record_ids = []
    for record in align:
        base = str(record.seq[col_idx])
        if base not in ["-", "."]:
            residue_count += 1
            if not base.isalpha() and base not in other_symbols:
                other_symbols.append(base)
                unusual_record_ids.append(record.id)
        else:
            gap_count += 1
    if unusual_record_ids:
        logger.debug(
            f"alignment column {col_idx} summary: residues={residue_count} gaps={gap_count} "
            f"unusual_symbols={other_symbols} unusual_record_ids={unusual_record_ids}"
        )
    return residue_count > 0


@cli.command("check-sprinx-headers")
@click.argument("fasta", type=click.Path(exists=True))
def check_sprinx_headers(fasta):
    """Fail fast if any header in FASTA can't be parsed for anticodon/aa by
    sprinx, rather than letting the whole sprinx run degrade silently."""
    from Bio import SeqIO
    from sprinx.common import header_to_aa, header_to_anticodon

    bad = []
    for record in SeqIO.parse(fasta, "fasta"):
        header = record.description
        if header_to_anticodon(header) is None or header_to_aa(header) is None:
            bad.append(record.id)

    if bad:
        raise ValueError(
            f"{len(bad)} record(s) in {fasta} have a header sprinx can't parse "
            f"for anticodon/aa (expected 'id|aa|anticodon|taxon', an "
            f"'anticodon=XXX' tag, or a GtRNAdb-style 'tRNA-{{AA}}-{{anticodon}}' "
            f"name): {', '.join(bad[:10])}"
            + (f", and {len(bad) - 10} more" if len(bad) > 10 else "")
        )


@cli.command("sprinx-to-seq-to-sprinzl")
@click.option("--output", required=True, help="Output FNAME (seq_to_sprinzl.tsv)")
@click.argument("sprinx_tsv", type=click.Path(exists=True))
def sprinx_to_seq_to_sprinzl(sprinx_tsv, output):
    """Convert sprinx's per-position sprinzl_mapping.tsv into the id/seq_pos/sprinzl
    format expected by ss_seq_to_sprinzl_final/ss_transform."""
    by_id = defaultdict(dict)   # seq_id -> {label: seq_pos}
    order = []                  # seq_id in order of first appearance
    labels = set()

    with open(sprinx_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            label = row["sprinzl_position"].strip()
            if not label:            # unlabeled position; nothing to map
                continue
            seq_id = row["seq_id"]
            if seq_id not in by_id:
                order.append(seq_id)
            # qutrna2 convention: 17A/20A/20B, not 17a/20a/20b; variable-arm
            # 'e' labels (e11, e21, ...) keep their lowercase 'e'.
            if not label.startswith(("e", "E")):
                label = label.upper()
            by_id[seq_id][label] = int(row["seq_index"]) + 1   # 1-indexed seq_pos
            labels.add(label)

    ordered_labels = sorted(labels, key=_label_sort_key)

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    with open(output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "seq_pos", "sprinzl"])
        for seq_id in order:
            positions = by_id[seq_id]
            for label in ordered_labels:
                w.writerow([seq_id, positions.get(label, "-"), label])

@cli.command("derive-sprinzl-labels")
@click.option("--output", required=True, help="Output FNAME (plain ordered label list)")
@click.argument("seq_to_sprinzl", type=click.Path(exists=True))
def derive_sprinzl_labels(seq_to_sprinzl, output):
    """Derive the flat, ordered Sprinzl label list plot_heatmap.R expects for
    --sprinzl from the union of labels present in a seq_to_sprinzl
    TSV, so modes that produce that TSV directly don't need a separately
    maintained labels file."""
    df = pd.read_csv(seq_to_sprinzl, sep="\t")
    labels = {label for label in df["sprinzl"].astype(str) if label not in ("-", ".")}
    ordered = sorted(labels, key=_label_sort_key)

    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    with open(output, "w") as f:
        for label in ordered:
            f.write(f"{label}\n")


@cli.command()
@click.option("--output", required=True, help="Output for aligned FASTA")
@click.argument("stk", type=click.Path(exists=True))
def stk_to_afasta(stk, output):
    align = AlignIO.read(stk, "stockholm")
    ss = str(align.column_annotations["secondary_structure"])

    id2seq = {}
    for record in align:
        new_seq = []
        for b, s in zip(str(record.seq), ss):
            if s == ".":
                if b != "-":
                    new_seq.append(b.lower())
                else:
                    new_seq.append("-")
            elif b == "-":
                new_seq.append(".")
            else:
                new_seq.append(b.upper())
        id2seq[record.id] = FastaRecord(record.id, "".join(new_seq))

    write_fasta_records(id2seq, output)


@cli.command()
@click.option("--consensus-labels", required=True, type=click.Path(exists=True))
@click.option("--output", required=True, help="Output FNAME")
@click.argument("afasta", type=click.Path(exists=True))
def afasta_to_sprinzl(afasta, consensus_labels, output):
    cl = pd.read_csv(consensus_labels, sep="\t")

    dfs = []
    faidx = pysam.FastaFile(afasta)
    for ref in faidx.references:
        seq = faidx[ref]
        la_sprinzl = []
        la_aln_pos = []
        la_seq_pos = []
        la_letter = []
        seq_pos = 0

        for aln_pos, (letter, label) in enumerate(zip(seq, cl["label"].to_list())):
            if letter in [".", "-"]:
                pass
            else:
                seq_pos += 1
                la_aln_pos.append(aln_pos)
                la_seq_pos.append(seq_pos)
                la_sprinzl.append(label)
                la_letter.append(letter)
        df = pd.DataFrame(
            {
                "id": ref,
                "seq_letter": la_letter,
                "seq_pos": la_seq_pos,
                "aln_pos": la_aln_pos,
                "sprinzl": la_sprinzl,
            }
        )
        dfs.append(df)

    pd.concat(dfs).to_csv(output, sep="\t", index=False, quoting=False)


@cli.command()
@click.option("--output", required=True, help="Output FNAME")
@click.option("--labels", required=True, type=click.Path(exists=True))
@click.argument("stk", type=click.Path(exists=True))
def consensus_labels(labels, stk, output):
    align = AlignIO.read(stk, "stockholm")
    ss = str(align.column_annotations["secondary_structure"])
    cl = pd.read_csv(labels, header=None)[0].to_list()

    aln_labels = []
    label_i = 0
    for s in ss:
        if s in ["(", ")", "<", ">", "{", "}"]:
            aln_labels.append(cl[label_i])
            label_i += 1
        elif s in [",", ":", "_", "-", "~"]:
            aln_labels.append(cl[label_i])
            label_i += 1
        elif s == ".":
            aln_labels.append("-")
        else:
            raise Exception(f"Unsupported secondary structure: {s}")
    try:
        assert len(aln_labels) == len(ss)
    except AssertionError:
        raise Exception("Mismatch of Sprinzl labels and secondary structure consensus alignment!")

    with open(output, "w") as f:
        f.write(f"aln_pos\tss\tlabel\n")
        for i, (s, l) in enumerate(zip(ss, aln_labels)):
            f.write(f"{i}\t{s}\t{l}\n")


@cli.command()
@click.option("--sprinzl", required=True, help="Sequence to Sprinzl.")
@click.option("--output", required=True, help="Output FNAME")
@click.option("--linker5", default=0, help="Length of 5' linker sequence")
@click.argument("jacusa2", type=click.Path(exists=True))
def transform(sprinzl, output, linker5, jacusa2):
    """Add Sprinzl coordinates to JACUS2A output"""

    sprinzl = pd.read_csv(sprinzl, sep="\t")
    jacusa = pd.read_csv(jacusa2, sep="\t")

    i = sprinzl["id"].isin(jacusa["trna"].unique())
    sprinzl = sprinzl.loc[i, ["id", "seq_pos", "sprinzl"]]
    sprinzl["seq_pos"] = sprinzl["seq_pos"].astype(str)

    jacusa["n_pos"] = jacusa["seq_position"] - linker5
    jacusa["n_pos"] = jacusa["n_pos"].astype(str)
    jacusa = (jacusa.merge(sprinzl,
                           how="left",
                           left_on=("trna", "n_pos" ),
                           right_on=("id", "seq_pos"),
                           indicator=True)
              .drop(columns=["n_pos", "seq_pos", "id"]))
    jacusa["sprinzl"] = jacusa["sprinzl"].fillna(".")

    jacusa.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    cli()
