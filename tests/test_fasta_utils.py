"""test_fasta_utils.py - unique_trna_seq and ignore_trnas on a small reference

input:   fasta_utils.Record objects built in the test
output:  pytest assertions
usage:   pytest tests/test_fasta_utils.py
env:     pysam and pandas, reached through the module under test
"""
import fasta_utils
import pytest

# five nuclear ids and one mitochondrial, the mix in a human reference. two
# pairs share a sequence and merge. one differs in the last id field. the
# other differs in the field before it
SEQUENCE_BY_ID = {"tRNA-Ala-AGC-1-1": "GGGGCCCC",
                  "tRNA-Ala-AGC-2-1": "GGGGCCCC",
                  "tRNA-Ala-AGC-3-1": "GGGGCCCA",
                  "tRNA-Gly-GCC-1-1": "TTTTAAAA",
                  "tRNA-Gly-GCC-1-2": "TTTTAAAA",
                  "MT-TC": "AAAATTTT"}


def reference():
    """A fresh id to record mapping, as it is built by read_records.

    Renaming is in place. Each test therefore builds its own.
    """
    return {name: fasta_utils.Record(name, seq)
            for name, seq in SEQUENCE_BY_ID.items()}


def test_one_record_survives_per_distinct_sequence():
    """Six ids, four distinct sequences, four records out. GGGGCCCA differs
    from GGGGCCCC by one base. It stays apart."""
    kept = fasta_utils.unique_trna_seq(reference())
    assert len(kept) == 4
    assert ({str(record.seq) for record in kept.values()}
            == set(SEQUENCE_BY_ID.values()))


def test_a_merged_record_takes_the_common_prefix_as_its_id():
    """tRNA-Ala-AGC-1-1 and -2-1 merge to tRNA-Ala-AGC. tRNA-Gly-GCC-1-1 and
    -1-2 merge to tRNA-Gly-GCC-1. The last character of the common prefix is
    dropped by normalize_trna_ids. That character is the separator where the
    ids differ in one field. The id of a unique record is kept."""
    assert set(fasta_utils.unique_trna_seq(reference())) == {
        "tRNA-Ala-AGC", "tRNA-Gly-GCC-1", "tRNA-Ala-AGC-3-1", "MT-TC"}


def test_ignore_removes_only_the_matching_trnas():
    """^MT is in three of the four data templates. The reference there is
    labelled by a euk CM. Those labels are for nuclear tRNAs alone."""
    kept = fasta_utils.ignore_trnas(reference(), "^MT")
    assert set(kept) == set(SEQUENCE_BY_ID) - {"MT-TC"}


def test_an_ignore_pattern_matching_every_id_raises():
    """^MT against a mito-only reference empties it. Under scheme mito the
    mt-tRNAs are labelled. The warning is in data_auto.yaml."""
    mito_only = {"MT-TC": fasta_utils.Record("MT-TC", "AAAATTTT")}
    with pytest.raises(fasta_utils.IgnorePatternError, match="No tRNAs left"):
        fasta_utils.ignore_trnas(mito_only, "^MT")


def test_an_ignore_pattern_matching_no_id_raises():
    """^mt in lowercase matches zero ids. re is case sensitive. Every
    mitochondrial tRNA would otherwise be kept by such a key."""
    with pytest.raises(fasta_utils.IgnorePatternError, match="Check pep file"):
        fasta_utils.ignore_trnas(reference(), "^mt")
