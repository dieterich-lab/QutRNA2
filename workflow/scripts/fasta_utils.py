import click
import enum
import gzip
import re
import pandas as pd
from collections import defaultdict
from os.path import commonprefix
from functools import partial

import pysam


class Record:
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq


class BaseChange(enum.Enum):
    U2T = enum.auto()
    T2U = enum.auto()


@click.group()
def cli():
    pass


def read_records(fname):
    faidx = pysam.Fastafile(fname)
    records = {ref: Record(ref, faidx[ref]) for ref in faidx.references}

    return records


def open_out(output):
    f = partial(gzip.open, mode="wt") if output.endswith("gz") else partial(open, mode="w")

    return f(output)


@cli.command()
@click.option("-o", "--output", required=True, type=click.Path())
@click.argument("FASTA", type=click.Path(exists=True))
def extract_seqids(fasta, output):
    """Extract tRNAs from fasta"""

    records = read_records(fasta)
    with open_out(output) as fout:
        for record_id in records:
            fout.write(f"{record_id}\n")


def do_remove_linker(record, linker5, linker3):
    if linker5 > 0:
        record.seq = record.seq[linker5:]

    if linker3 > 0:
        record.seq = record.seq[:-linker3]

    return record


def do_base_change(record, base1, base2):
    record.seq = record.seq.replace(base1, base2)

    return record


def do_reverse(record, linker5=0, linker3=0):
    record.id = f"{record.id}_rev"
    s = record.seq

    start = linker5
    end = len(s) - linker3
    l5 = ""
    if linker5 > 0:
        l5 = s[0:start]
    l3 = ""
    if linker3 > 0:
        l3 = s[(len(s) - linker3):len(s)]
    trna = s[start:end][::-1]
    record.seq = f"{l5}{trna}{l3}"

    return record


@cli.command()
@click.option("-5", "--linker5", type=str)
@click.option("-3", "--linker3", type=str)
@click.option("--linker5-length", type=int, default=0)
@click.option("--linker3-length", type=int, default=0)
@click.option("--remove-linker5", type=int)
@click.option("--remove-linker3", type=int)
@click.option("-c", "--base-change", type = click.Choice(BaseChange, case_sensitive=False), help="Base change: U2T or T2U.")
@click.option("-i", "--ignore", type=str)
@click.option("-r", "--reverse", is_flag=True, default=False, help="Reverse sequence.")
@click.option("-u", "--unique-seq", is_flag=True, default=False, help="Only unique sequences.")
@click.option("-o", "--output", required=True, type=click.Path())
@click.argument("FASTA", type=click.Path())
def transform(fasta, linker5, linker3, linker5_length, linker3_length, remove_linker5, remove_linker3, base_change, ignore, reverse, unique_seq, output):
    """Transform fasta"""

    tasks = []
    if remove_linker5 or remove_linker3:
        def helper(record):
            return do_remove_linker(record, remove_linker5, remove_linker3)
        tasks.append(helper)

    if base_change:
        bc = str(base_change).split("2")
        def helper(record):
            return do_base_change(record, bc[0], bc[1])
        tasks.append(helper)

    if reverse:
        def helper(record):
            return do_reverse(record, linker5_length, linker3_length)
        tasks.append(helper)


    if linker5:
      def helper(record):
        record.seq = f"{linker5}{str(record.seq)}"

        return(record)
      tasks.append(helper)
    if linker3:
      def helper(record):
        record.seq = f"{record.seq}{linker3}"

        return(record)
      tasks.append(helper)

    records = read_records(fasta)
    if ignore:
      records = ignore_trnas(records, ignore)

    if unique_seq:
      records = unique_trna_seq(records)

    with open_out(output) as fout:
        for record in records.values():
            for task in tasks:
                record = task(record)
            fout.write(f">{record.id}\n")
            fout.write(f"{record.seq}\n")


def ignore_trnas(records, ignore):
    """ Remove trnas from fasta"""

    filtered = dict(records)

    pattern = re.compile(rf"{ignore}")
    trnas = [trna for trna in filtered if pattern.search(trna)]

    # remove trnas such as defined by ignore
    for trna in trnas:
        del filtered[trna]

    if not filtered:
        raise Exception("No tRNAs left!")

    if ignore and len(filtered) == len(records):
        raise Exception("No tRNAs were filtered. Check pep file!")

    return filtered


def normalize_trna_ids(trna_ids: list):
  if len(trna_ids) == 1:
    return trna_ids[0]

  common = commonprefix(trna_ids)

  # add suffix FIXME

  return common[:-1]


def unique_trna_seq(records):
  trna_seq_to_id = defaultdict(set)
  for trna_id, record in records.items():
    trna_seq_to_id[str(record.seq)].add(trna_id)

  # normalize labels
  new_records = {}
  for trna_seq, trna_ids in trna_seq_to_id.items():
    trna_id = normalize_trna_ids(list(trna_ids))
    record = records[list(trna_ids)[0]]
    record.id = trna_id
    record.name = trna_id
    record.description = ""
    new_records[trna_id] = record

  return new_records


@cli.command()
@click.option("-o", "--output", type=click.Path())
@click.argument("FASTA", type=click.Path(exists=True))
def infer_annotation(fasta, output):
  """Infer tRNA annotation from fasta"""

  trnas = read_records(fasta)

  df = pd.DataFrame.from_dict({"trna": trnas.keys()})

  df["type"] = "unknown"
  is_mt = df["trna"].str.startswith("MT")
  df["type"][is_mt] = "mt"
  df["type"][~is_mt] = "nuclear"

  if df["trna"].str.contains("tRNA").any():
    info = df["trna"].str.extract(r'.*tRNA-(?P<amino_acid>[^-]+)-(?P<anti_codon>[A-Z]+)')
  else:
    info = df["trna"].str.extract(r'(?P<amino_acid>[^-]+)-(?P<anti_codon>[A-Z]+)')

  df = pd.concat([df, info], axis=1)
  df["seq"] = [str(record.seq).upper().replace("T", "U") for record in trnas.values()]

  df.to_csv(output, index=False, sep="\t")


if __name__ == "__main__":
    cli()
