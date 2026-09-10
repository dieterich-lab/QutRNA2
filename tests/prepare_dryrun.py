"""Build a working directory for a snakemake dry run.

Copies the shipped templates from examples/ so a broken one fails here, and
synthesises only what a dry run cannot supply: a reference FASTA and empty BAMs.

usage: python tests/prepare_dryrun.py <workdir>
"""
import csv
import os
import shutil
import sys

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# human mt-tRNAs, from sprinx data/mito/canonical.fa (mtdbD00063498, D00063499)
REFERENCE = {
    "Homo_sapiens_mt-tRNA-Phe-GAA-1-1":
        "GTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACA",
    "Homo_sapiens_mt-tRNA-Val-TAC-1-1":
        "CAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGA",
}

# RNA004 adapters, lengths matching linker5/linker3 in the data template.
# remove_linker strips these before sprinx sees the tRNA.
LINKER5 = "CCTAAGAGCAAGAAGAAGCCTGGN"
LINKER3 = "GGCTTCTTCTTGCTCTTAGGAAAAAAAAAAAA"


def _disable_random_alignment(text):
    """that filter reads mapping stats, which existing_mapping never produces."""
    out, in_block = [], False
    for line in text.splitlines(keepends=True):
        if line.startswith("  random_alignment:"):
            in_block = True
        elif in_block and line.strip().startswith("apply:"):
            line, in_block = line.replace("True", "False"), False
        out.append(line)
    return "".join(out)


def main(workdir):
    data = os.path.join(workdir, "data")
    os.makedirs(data, exist_ok=True)

    shutil.copy(os.path.join(REPO, "examples", "data", "data_auto.yaml"),
                os.path.join(data, "data.yaml"))
    shutil.copy(os.path.join(REPO, "examples", "sample_table_bam.tsv"),
                os.path.join(data, "sample_table.tsv"))

    analysis = open(os.path.join(REPO, "examples", "analysis",
                                 "existing_mapping.yaml")).read()
    with open(os.path.join(workdir, "analysis.yaml"), "w") as fh:
        fh.write(_disable_random_alignment(analysis))

    with open(os.path.join(workdir, "reference.fasta"), "w") as fh:
        for header, seq in REFERENCE.items():
            fh.write(f">{header}\n{LINKER5}{seq}{LINKER3}\n")

    with open(os.path.join(data, "sample_table.tsv"), newline="") as fh:
        bams = {row["bam"] for row in csv.DictReader(fh, delimiter="\t")}
    for bam in bams:
        open(os.path.join(workdir, bam), "w").close()

    print(f"dry-run workdir ready: {workdir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")