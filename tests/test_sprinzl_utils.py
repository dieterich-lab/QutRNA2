"""test_sprinzl_utils.py - label ordering and the sprinx TSV conversion

input:   TSV fixtures written under tmp_path
output:  pytest assertions
usage:   pytest tests/test_sprinzl_utils.py
env:     pysam and pandas, plus biopython for the Stockholm reader
"""
import csv

import pytest
import sprinzl_utils


def test_labels_sort_in_5p_to_3p_order():
    """One pass over every label shape the pipeline emits. The awkward cases:
    60Z precedes 60AA, the variable arm falls inside the 45/46 gap, and its 3'
    stem counts down (e27 before e21) because that is how it reads 5'->3'."""
    shuffled = ["46", "e21", "60AA", "17A", "9", "e1", "60Z", "20B", "45",
                "17", "e27", "60", "76", "e11", "20", "e5", "60A", "20A", "e17"]
    assert sorted(shuffled, key=sprinzl_utils._label_sort_key) == [
        "9", "17", "17A", "20", "20A", "20B", "45",
        "e11", "e17", "e1", "e5", "e27", "e21",
        "46", "60", "60A", "60Z", "60AA", "76"]


def test_an_unrecognised_label_is_named_in_the_error():
    """_label_sort_key('ins3') raises ValueError. 'ins3' is in the message.
    Label shapes keep arriving from sprinx."""
    with pytest.raises(ValueError, match="ins3"):
        sprinzl_utils._label_sort_key("ins3")


def test_sprinx_conversion_normalizes_case_and_indexing(tmp_path):
    """0-based indices and lowercase insertion codes are written by sprinx.
    1-based seq_pos and 17A/20B are expected downstream. An 'e' label is left
    alone. Only a labelled position gets a row."""
    src, out = tmp_path / "sprinx.tsv", tmp_path / "seq_to_sprinzl.tsv"
    with open(src, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["seq_id", "seq_index", "nucleotide", "sprinzl_position"])
        w.writerows([["tRNA-Ala", 0, "G", "17"],
                     ["tRNA-Ala", 1, "A", "17a"],
                     ["tRNA-Ala", 2, "C", "20b"],
                     ["tRNA-Ala", 3, "U", "e11"],
                     ["tRNA-Ala", 4, "G", ""]])

    sprinzl_utils.sprinx_to_seq_to_sprinzl.callback(str(src), str(out))

    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["sprinzl"] for r in rows} == {"17", "17A", "20B", "e11"}
    assert next(r["seq_pos"] for r in rows if r["sprinzl"] == "17") == "1"


def test_derived_labels_are_ordered_and_skip_blanks(tmp_path):
    """A '-' blanks a position out of the plot."""
    src, out = tmp_path / "seq_to_sprinzl.tsv", tmp_path / "labels.txt"
    with open(src, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "seq_pos", "sprinzl"])
        w.writerows([["tRNA-Ala", 1, "20A"], ["tRNA-Ala", 2, "17"],
                     ["tRNA-Ala", 3, "-"], ["tRNA-Ala", 4, "20"]])

    sprinzl_utils.derive_sprinzl_labels.callback(str(src), str(out))

    with open(out) as fh:
        assert fh.read().split() == ["17", "20", "20A"]
