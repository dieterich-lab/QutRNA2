import click
import pandas as pd
import pysam
from Bio import AlignIO


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
      {"id": ref,
      "seq_letter": la_letter,
      "seq_pos": la_seq_pos,
      "aln_pos": la_aln_pos,
      "sprinzl": la_sprinzl})
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
