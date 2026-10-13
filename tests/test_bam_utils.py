"""test_bam_utils.py - trim_cigar_helper and do_trim_cigar on synthetic cigars

input:   cigar tuples, and pysam records built without a header
output:  pytest assertions
usage:   pytest tests/test_bam_utils.py
env:     pysam, reached through the module under test
"""
import bam_utils
import pysam
import pytest

# SAM spec section 1.4, the "consumes query" and "consumes reference" columns
CONSUMES_REFERENCE = {pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP,
                      pysam.CEQUAL, pysam.CDIFF}
CONSUMES_QUERY = {pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP,
                  pysam.CEQUAL, pysam.CDIFF}

# the walk stops at an aligned or clipped operation
REWRITTEN_OPERATIONS = [pysam.CDIFF, pysam.CDEL, pysam.CINS]


def read(cigartuples, reference_start=0):
    """A pysam record carrying `cigartuples`, built without a header.

    cigartuples is read back as a list. The comparison in do_trim_cigar is
    against that list.
    """
    record = pysam.AlignedSegment()
    record.query_name = "read"
    record.reference_start = reference_start
    record.cigartuples = cigartuples
    return record


def cigar_starting_with(operation, length=3):
    """A fresh cigar of one leading `operation` of `length`, then 47M.

    The 47M is an arbitrary aligned run.
    """
    return [(operation, length), (pysam.CMATCH, 47)]


def bases_consumed(cigar, operations):
    """The base count for `cigar`, restricted to `operations`."""
    # an empty tuple is left in the cigar by a dropped deletion
    return sum(length for op, length in filter(None, cigar) if op in operations)


@pytest.mark.parametrize("operation", REWRITTEN_OPERATIONS)
def test_the_offset_counts_the_reference_bases_removed(operation):
    """3X47M and 3D47M give offset 3. 3I47M gives 0. Query alone is consumed
    by I."""
    _, offset = bam_utils.trim_cigar_helper(cigar_starting_with(operation, 3))
    assert offset == (3 if operation in CONSUMES_REFERENCE else 0)


@pytest.mark.parametrize("operation", REWRITTEN_OPERATIONS)
def test_trimming_preserves_the_query_length(operation):
    """3I47M becomes 3S47M. 3D47M becomes 47M."""
    trimmed, _ = bam_utils.trim_cigar_helper(cigar_starting_with(operation))
    assert (bases_consumed(trimmed, CONSUMES_QUERY)
            == bases_consumed(cigar_starting_with(operation), CONSUMES_QUERY))


@pytest.mark.parametrize("operation", [pysam.CMATCH, pysam.CEQUAL,
                                       pysam.CSOFT_CLIP, pysam.CHARD_CLIP])
def test_a_deletion_inside_the_alignment_survives(operation):
    """The 3D in 5M3D47M is kept. Only a leading indel is a mapper artifact."""
    first_operation = [(operation, 5)]
    trimmed, offset = bam_utils.trim_cigar_helper(
        first_operation + cigar_starting_with(pysam.CDEL))
    assert trimmed == first_operation + cigar_starting_with(pysam.CDEL)
    assert offset == 0


def test_an_untouched_read_is_reported_as_untouched():
    """2S10M is already clean. False is returned and reference_start is left
    alone. A wrong comparison would count it as trimmed."""
    record = read([(pysam.CSOFT_CLIP, 2), (pysam.CMATCH, 10)])
    assert bam_utils.do_trim_cigar(record) is False
    assert record.reference_start == 0


def test_the_offset_reaches_the_read():
    """2X3S10M trims to 5S10M. The 2X becomes a soft clip and merges with the
    3S. reference_start advances by 2."""
    record = read([(pysam.CDIFF, 2), (pysam.CSOFT_CLIP, 3),
                   (pysam.CMATCH, 10)], reference_start=100)
    assert bam_utils.do_trim_cigar(record) is True
    assert record.reference_start == 102
    assert record.cigartuples == [(pysam.CSOFT_CLIP, 5), (pysam.CMATCH, 10)]
