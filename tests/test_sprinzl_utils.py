"""Unit tests for workflow/scripts/sprinzl_utils.py.

Only the parts needing no external tool: label ordering, which decides heatmap
column order, and the sprinx TSV conversion feeding it.
"""
import csv
import importlib.util
import os
import sys

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")


def _load(name):
    path = os.path.join(SCRIPTS, f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sprinzl_utils = _load("sprinzl_utils")


def test_labels_sort_in_5p_to_3p_order():
    """One pass over every label shape the pipeline emits. The awkward cases:
    60Z precedes 60AA, the variable arm sits inside the 45/46 gap, and its 3'
    stem counts down (e27 before e21) because that is how it reads 5'->3'."""
    shuffled = ["46", "e21", "60AA", "17A", "9", "e1", "60Z", "20B", "45",
                "17", "e27", "60", "76", "e11", "20", "e5", "60A", "20A", "e17"]
    assert sorted(shuffled, key=sprinzl_utils._label_sort_key) == [
        "9", "17", "17A", "20", "20A", "20B", "45",
        "e11", "e17", "e1", "e5", "e27", "e21",
        "46", "60", "60A", "60Z", "60AA", "76"]


def test_sprinx_conversion_normalizes_case_and_indexing(tmp_path):
    """sprinx writes 0-based indices and lowercase insertion codes; downstream
    wants 1-based seq_pos and 17A/20B, with 'e' labels left alone. Positions
    sprinx could not label carry no row."""
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

    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert {r["sprinzl"] for r in rows} == {"17", "17A", "20B", "e11"}
    assert next(r["seq_pos"] for r in rows if r["sprinzl"] == "17") == "1"


def test_derived_labels_are_ordered_and_skip_blanks(tmp_path):
    """A '-' blanks a position out of the plot rather than naming a coordinate."""
    src, out = tmp_path / "seq_to_sprinzl.tsv", tmp_path / "labels.txt"
    with open(src, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "seq_pos", "sprinzl"])
        w.writerows([["tRNA-Ala", 1, "20A"], ["tRNA-Ala", 2, "17"],
                     ["tRNA-Ala", 3, "-"], ["tRNA-Ala", 4, "20"]])

    sprinzl_utils.derive_sprinzl_labels.callback(str(src), str(out))

    assert open(out).read().split() == ["17", "20", "20A"]