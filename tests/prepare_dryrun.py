"""prepare_dryrun.py - build a working directory for one snakemake dry run

input:   one data template and one analysis template from examples/
output:  WORKDIR holding analysis.yaml, data/, and the synthesised inputs
usage:   python tests/prepare_dryrun.py WORKDIR [--data NAME] [--analysis NAME]
env:     pyyaml
notes:   templates are copied verbatim. A broken one fails here.
         Synthesised: the reference FASTA, empty read files, Sprinzl inputs
"""
import argparse
import csv
import os
import pathlib
import shutil

import yaml

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def templates(kind):
    """Stems of the templates shipped under examples/<kind>/."""
    listing = os.listdir(os.path.join(REPO, "examples", kind))
    return sorted(name[:-5] for name in listing if name.endswith(".yaml"))


def combinations():
    """Every data x analysis pair, then minimal once per analysis template.

    The data templates differ only in how they supply Sprinzl coordinates.
    minimal unsets coords, which drops the Sprinzl rules from the DAG. One
    data template therefore covers the whole minimal arm.
    """
    for data in templates("data"):
        for analysis in templates("analysis"):
            yield data, analysis, ""
    for analysis in templates("analysis"):
        yield templates("data")[0], analysis, "--minimal"


# human mt-tRNAs, from sprinx data/mito/canonical.fa (mtdbD00063498, D00063499)
REFERENCE = {
    "Homo_sapiens_mt-tRNA-Phe-GAA-1-1":
        "GTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACA",
    "Homo_sapiens_mt-tRNA-Val-TAC-1-1":
        "CAGAGTGTAGCTTAACACAAAGCACCCAACTTACACTTAGGAGATTTCAACTTAACTTGACCGCTCTGA",
}

# RNA004 adapters, lengths matching linker5/linker3 in the data templates.
# remove_linker strips these before sprinx sees the tRNA
LINKER5 = "CCTAAGAGCAAGAAGAAGCCTGGN"
LINKER3 = "GGCTTCTTCTTGCTCTTAGGAAAAAAAAAAAA"

# random_alignment reads mapping stats. These runs do not produce them
NEEDS_NO_RANDOM_ALIGNMENT = ("existing_mapping",)


def write_reference(workdir):
    path = os.path.join(workdir, "reference.fasta")
    with open(path, "w") as fh:
        fh.writelines(f">{header}\n{LINKER5}{seq}{LINKER3}\n"
                      for header, seq in REFERENCE.items())


def write_afasta(workdir, pep):
    """An aligned FASTA plus one consensus label per column.

    The template ships both keys empty. Padding to equal length is enough.
    """
    width = max(len(s) for s in REFERENCE.values())
    with open(os.path.join(workdir, "afasta.fasta"), "w") as fh:
        fh.writelines(f">{header}\n{seq}{'-' * (width - len(seq))}\n"
                      for header, seq in REFERENCE.items())
    with open(os.path.join(workdir, "consensus_labels.txt"), "w") as fh:
        fh.write("\n".join(str(i) for i in range(1, width + 1)) + "\n")

    pep["qutrna2"]["sprinzl"]["afasta"] = "afasta.fasta"
    pep["qutrna2"]["sprinzl"]["consensus_labels"] = "consensus_labels.txt"


def prepare_data(workdir, name, *, minimal=False):
    """Copy the data template and supply the files named by its sprinzl block.

    minimal drops coords, which every shipped template sets.
    """
    with open(os.path.join(REPO, "examples", "data", f"{name}.yaml")) as fh:
        pep = yaml.safe_load(fh)
    if minimal:
        pep["qutrna2"].pop("coords", None)
    sprinzl = pep["qutrna2"].get("sprinzl", {})

    # the pep names these by path; the curated copies are in the repo's data/
    for key in ("labels", "cm", "seq_to_sprinzl"):
        if sprinzl.get(key):
            target = os.path.join(workdir, sprinzl[key])
            os.makedirs(os.path.dirname(target) or workdir, exist_ok=True)
            shutil.copy(os.path.join(REPO, "data",
                                     os.path.basename(sprinzl[key])), target)
    if name == "data_afasta":
        write_afasta(workdir, pep)

    with open(os.path.join(workdir, "data", "data.yaml"), "w") as fh:
        yaml.safe_dump(pep, fh, sort_keys=False)


def prepare_analysis(workdir, name, *, minimal=False):
    """minimal drops plots, which every shipped template sets."""
    with open(os.path.join(REPO, "examples", "analysis", f"{name}.yaml")) as fh:
        config = yaml.safe_load(fh)
    if name in NEEDS_NO_RANDOM_ALIGNMENT:
        config["filter"]["random_alignment"]["apply"] = False
    if minimal:
        config.pop("plots", None)
    with open(os.path.join(workdir, "analysis.yaml"), "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)


def prepare_reads(workdir, analysis):
    """Copy the sample table implied by the analysis template, then touch
    its files."""
    column = "bam" if analysis in NEEDS_NO_RANDOM_ALIGNMENT else "fastq"
    table = os.path.join(workdir, "data", "sample_table.tsv")
    shutil.copy(os.path.join(REPO, "examples", f"sample_table_{column}.tsv"), table)

    with open(table, newline="") as fh:
        names = {row[column] for row in csv.DictReader(fh, delimiter="\t")}
    for name in names:
        pathlib.Path(workdir, name).touch()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("workdir", nargs="?")
    parser.add_argument("--data", choices=templates("data"), default="data_auto")
    parser.add_argument("--analysis", choices=templates("analysis"),
                        default="existing_mapping")
    parser.add_argument("--minimal", action="store_true",
                        help="drop coords and plots, leaving the defaults to supply them")
    parser.add_argument("--list", action="store_true",
                        help="print every combination, one per line, and exit")
    args = parser.parse_args()

    if args.list:
        for data, analysis, extra in combinations():
            print(data, analysis, extra)
        return
    if not args.workdir:
        parser.error("workdir is required")

    os.makedirs(os.path.join(args.workdir, "data"), exist_ok=True)
    prepare_data(args.workdir, args.data, minimal=args.minimal)
    prepare_analysis(args.workdir, args.analysis, minimal=args.minimal)
    prepare_reads(args.workdir, args.analysis)
    write_reference(args.workdir)

    print(f"dry-run workdir ready: {args.workdir} ({args.data} + {args.analysis})")


if __name__ == "__main__":
    main()
