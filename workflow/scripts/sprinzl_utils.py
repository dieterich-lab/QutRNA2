import os
import logging
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


# FIXME

# Organism/model-specific Sprinzl schemes (from tDRnamer tdrdbutils.py).
# Labels are used in order for non-gap consensus structure columns.
SPRINZL_SCHEMES = {
    "euk": [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
        "17", "17a", "18", "19", "20", "20a", "20b", "21", "22", "23", "24", "25", "26", "27", "28", "29",
        "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45",
        "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e1", "e2", "e3", "e4", "e5",
        "e27", "e26", "e25", "e24", "e23", "e22", "e21",
        "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
        "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76"
    ],
    "arch": [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
        "17", "17a", "18", "19", "20", "20a", "20b", "21", "22", "23", "24", "25", "26", "27", "28", "29",
        "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45",
        "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e1", "e2", "e3", "e4",
        "e27", "e26", "e25", "e24", "e23", "e22", "e21",
        "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
        "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76"
    ],
    "bact": [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
        "17", "17a", "18", "19", "20", "20a", "20b", "21", "22", "23", "24", "25", "26", "27", "28", "29",
        "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45",
        "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e1", "e2", "e3", "e4",
        "e27", "e26", "e25", "e24", "e23", "e22", "e21",
        "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
        "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76"
    ],
    "mito": [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
        "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37",
        "38", "39", "40", "41", "42", "43", "44", "45",
        "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e1", "e2", "e3", "e4",
        "e27", "e26", "e25", "e24", "e23", "e22", "e21",
        "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64",
        "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76"
    ]
}


def _is_label_column(ss_char):
    return ss_char != "."


def _normalize_sprinzl_label(label):
    s = str(label)
    if s.startswith("e"):
        return s
    if len(s) > 1 and s[-1].isalpha():
        return f"{s[:-1]}{s[-1].upper()}"
    return s


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


def _available_label_columns(align, ss):
    return [
        i
        for i, ss_char in enumerate(ss)
        if _is_label_column(ss_char) and _column_has_residue(align, i)
    ]


def _labels_from_scheme_and_alignment(align, ss, scheme_labels):
    label_columns = [i for i, ss_char in enumerate(ss) if _is_label_column(ss_char)]
    if len(label_columns) > len(scheme_labels):
        raise ValueError(
            f"Scheme has {len(scheme_labels)} labels but alignment needs {len(label_columns)} label columns."
        )

    available_label_columns = set(_available_label_columns(align, ss))
    labels = []
    for label_idx, col_idx in enumerate(label_columns):
        label = _normalize_sprinzl_label(scheme_labels[label_idx])
        if _column_has_residue(align, col_idx):
            labels.append(label)
        else:
            labels.append("-")

    labeled_available_columns = [
        col_idx
        for col_idx, label in zip(label_columns, labels)
        if col_idx in available_label_columns and label != "-"
    ]
    if len(labeled_available_columns) != len(available_label_columns):
        missing_columns = sorted(available_label_columns.difference(labeled_available_columns))
        raise ValueError(
            "Auto-label generation failed: some biologically available label columns were not assigned "
            f"a Sprinzl label ({missing_columns})."
        )

    return labels

@cli.command()
@click.option("--output", required=True, help="Output labels file.")
@click.option(
    "--scheme",
    required=True,
    type=click.Choice(["euk", "arch", "bact", "mito"], case_sensitive=False),
    help="Sprinzl scheme: euk|arch|bact|mito",
)
@click.option("--debug", is_flag=True, help="Enable debug logging for residue checks.")
@click.argument("stk", type=click.Path(exists=True))
def auto_labels(stk, output, scheme, debug):
    """Generate Sprinzl labels from CM alignment consensus structure."""
    if debug:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s:%(name)s:%(message)s")

    try:
        align = AlignIO.read(stk, "stockholm")
    except Exception as e:
        raise IOError(f"Failed to read alignment from {stk}: {e}") from e

    if "secondary_structure" not in align.column_annotations:
        raise ValueError(f"{stk}: missing secondary_structure annotation (invalid cmalign output?)")

    ss = str(align.column_annotations["secondary_structure"])
    selected_scheme = scheme.lower()

    if selected_scheme not in SPRINZL_SCHEMES:
        raise ValueError(
            f"Invalid Sprinzl scheme '{scheme}'. "
            "Specify one of: euk, arch, bact, mito."
        )
    labels = _labels_from_scheme_and_alignment(align, ss, SPRINZL_SCHEMES[selected_scheme])
    required = sum(1 for c in ss if _is_label_column(c))
    if len(labels) != required:
        raise ValueError(
        f"Auto-label generation failed: generated {len(labels)} labels "
        f"for {required} label columns."
        )

    # Write labels
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
    try:
        with open(output, "w", encoding="utf-8") as f:
            for label in labels:
                f.write(f"{label}\n")
    except Exception as e:
        raise IOError(f"Failed to write {output}: {e}") from e

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
