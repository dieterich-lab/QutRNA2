import os
import shlex

import pandas as pd
from snakemake.io import directory, unpack

global DEFAULT_SCORE
global REF_FASTA
global READS
global SCORES
global SPRINZL_LABELS
global TRNA_ANNOTATION
global MAX_SCORES


################################################################################
# Plot JACUSA2 score: Mis+Del+Ins
################################################################################


def get_read_feature_plots(ftype):
  # return container
  plots = {}

  features = ["record_count", "read_length"]
  types = ["condition", "sample", "subsample"]
  for feature in features:
    plots.setdefault(feature, {})
    # default overview plots of read features
    for type in types:
      plots[feature][type] = f"results/plots/{feature}/{type}.{ftype}"
    # custom read feature specific plots

    plots[feature]["custom"] = {}
    for plot in config["plots"].get(feature, []):
      plots[feature]["custom"][plot["id"]] = f"results/plots/{feature}/custom/{plot['id']}.{ftype}"

  if READS == "fastq":
    plots["threshold_summary"] = {
      "default": f"results/plots/alignment/threshold_summary.{ftype}"
    }
    plots["threshold_summary"]["custom"] = {}
    for plot in config["plots"].get("threshold_summary", []):
      plots["threshold_summary"]["custom"] = f"results/plots/alignment/threshold_summary/custom/{plot['id']}.{ftype}"

  # remove keys with empty values
  plots = {k: v for k, v in plots.items() if v}

  return plots


def get_heatmap_plots(ftype):
  plots = {}

  contrasts = pep.config["qutrna2"]["contrasts"]
  bam_types = ["final"]
  if config["call_filtered"]:
    for f in FILTERS_APPLIED:
      bam_types.append(f"filtered-{f}")
  for contrast in contrasts:
    cond1 = contrast["cond1"]
    cond2 = contrast["cond2"]
    for plot in subsampled_plots(cond1, cond2):
      plot_id = plot["id"]
      for bam_type in bam_types:
        fname = f"results/plots/scores/cond1~{cond1}/cond2~{cond2}/{plot_id}/bam~{bam_type}/files.tsv"
        df = pd.read_csv(fname, sep="\t")
        df = df[df["type"] == ftype]
        for row in df.itertuples(index=False):
          # the stem already carries the group as a prefix; keying by group
          # alone collapses heatmap, heatmap_cov and heatmap_mods onto one entry
          stem = os.path.splitext(os.path.basename(row.name))[0]
          plots.setdefault(cond1, {}).setdefault(cond2, {}).setdefault(plot_id, {}).setdefault(bam_type, {})[stem] = row.name

  return plots


def get_plot(plot_type, plot_id):
  plots = config["plots"][plot_type]
  for plot in plots:
    if plot["id"] == plot_id:
      return plot

  return {}


def _plot_heatmap_opts(wildcards, input):
  plot = get_plot("heatmap", wildcards.plot_id)

  opts = {}

  col_mapping = {
      "seq": "seq_position",
      "sprinzl": "sprinzl"}
  opts["--position_column"] = col_mapping[pep.config["qutrna2"]["coords"]]

  if "group" in plot:
    opts[f"--group"] =  plot["group"]
  if plot.get("title"):
    opts["--title"] = f"'{plot['title']}'"

  linker5 = str(pep.config["qutrna2"]["linker5"])
  opts["--five_adapter"] = linker5

  for contrast in pep.config["qutrna2"]["contrasts"]:
    if wildcards.COND1 == contrast["cond1"] and wildcards.COND2 == contrast["cond2"] and "flag" in contrast:
      flag_positions = ",".join(contrast["flag"])
      opts["--flag_positions"] =  f"'{flag_positions}'"

  score = DEFAULT_SCORE
  if "score" in plot:
    score = plot.get("score", DEFAULT_SCORE)
  opts["--score_column"] = score.split("::")[0]
  if "sprinzl" in pep.config["qutrna2"]:
    opts["--sprinzl"] = input.sprinzl

  if "abbrevs" in pep.config["qutrna2"].get("mods", {}):
    opts["--abbrevs"] = input.abbrevs

  if plot.get("harmonize_scaling"):
    if plot["harmonize_scaling"] in ["CONTRASTS", "OVERALL"]:
      opts["--max_scores"] = input.max_scores
    opts["--harmonize_scaling"] = plot["harmonize_scaling"]

  s = " ".join([f"{key}={value}" for key, value in opts.items()])
  if "opts" in plot and plot["opts"]:
    s = f"{s} {plot['opts']}"

  return s


def _plot_heatmap_input(wildcards):
  plot = get_plot("heatmap", wildcards.plot_id)
  input = {
    "scores": f"results/jacusa2/cond1~{{COND1}}/cond2~{{COND2}}/bam~{{bam_type}}/{SCORES}",
    "trna_annotation": TRNA_ANNOTATION
  }

  if "sprinzl" in pep.config["qutrna2"]:
    input["sprinzl"] = SPRINZL_LABELS

  if "abbrevs" in pep.config["qutrna2"].get("mods",{}):
    global MOD_ABBREVS

    input["abbrevs"] = MOD_ABBREVS

  if plot.get("harmonize_scaling", "") in ["CONTRASTS", "OVERALL"]:
    input["max_scores"] = MAX_SCORES

  return input


checkpoint plot_heatmap:
  input: unpack(_plot_heatmap_input)
  output: dir=directory("results/plots/scores/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}"),
          files="results/plots/scores/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}/files.tsv"
  conda: "qutrna2"
  log: "logs/plot/heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}.log"
  params: opts=_plot_heatmap_opts,
          basedir=workflow.basedir
  shell: """
    ( export QUTRNA2="{params.basedir}"
      mkdir -p {output.dir} && \
      Rscript {workflow.basedir}/scripts/plot_heatmap.R \
          --condition1 '{wildcards.COND1}' \
          --condition2 '{wildcards.COND2}' \
          --bam_type '{wildcards.bam_type}' \
          --trna_annotation {input.trna_annotation:q} \
          --output_dir {output.dir:q} \
          {params.opts} \
          {input.scores:q} \
          ) 2> {log:q}
  """


def _plot_normalized_heatmap_input(wildcards):
  input = {"perm": "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/permutation.tsv"}
  if COORDS == "sprinzl":
    input["sprinzl"] = "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/scores_sprinzl.tsv"

  return input


rule plot_normalized_heatmap:
  input: unpack(_plot_normalized_heatmap_input)
  output: pdf="results/plots/normalized/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/normalized_heatmap.pdf",
          png="results/plots/normalized/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/normalized_heatmap.png"
  conda: "qutrna2"
  log: "logs/plot/normalized_heatmap/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}.log"
  params: sprinzl=lambda wildcards, input: f"--sprinzl {input.sprinzl}" if "sprinzl" in input.keys() else "",
          fdr=config["jacusa2"]["fdr"],
          opts=lambda wildcards: _norm_opts(config["plots"].get("all_trnas", {})),
          linker5=pep.config["qutrna2"]["linker5"]
  shell: """
    python {workflow.basedir:q}/scripts/plot_norm_heatmap.py \
        --conditions '{wildcards.COND1},{wildcards.COND2}' \
        --fdr {params.fdr} \
        --five_adapter {params.linker5} \
        {params.opts} \
        --out {output.pdf:q} \
        {params.sprinzl} \
        {input.perm:q} > {log:q} 2>&1
  """


def _norm_opts(plot):
  # the entry's plot_heatmap.R flags, under plot_norm_heatmap.py's names
  flags = {"--trnas": "--trnas", "--mark_positions": "--mark",
           "--min_reads": "--min_reads", "--remove_prefix": "--remove_prefix"}
  opts = [f"{new}={shlex.quote(value)}" for old, new in flags.items()
          if (value := heatmap_opt(plot, old))]
  # a bare flag, which heatmap_opt does not read
  if "--hide_varm" in shlex.split(plot.get("opts") or ""):
    opts.append("--hide_varm")
  if plot.get("coverage_split"):
    opts.append("--coverage_split")
  elif plot.get("group"):
    opts.append(f"--group={plot['group']} --trna_annotation={TRNA_ANNOTATION}")
  scaling = str(plot.get("harmonize_scaling", "NONE"))
  if scaling.replace(".", "", 1).isdigit():
    opts.append(f"--max_score={scaling}")
  return " ".join([f"--score_column={plot.get('norm_score', DEFAULT_NORM_SCORE)}", *opts])


def _normalized_set_opts(wildcards):
  opts = _norm_opts(get_plot("heatmap", wildcards.plot_id))
  if (tables := _scale_tables(wildcards)):
    opts += " " + " ".join(["--scale_with", *tables])
  return opts


def _scale_tables(wildcards):
  """Permutation tables over which the colour scale is shared: the entry's
  contrasts for harmonize_scaling CONTRASTS, every contrast for OVERALL."""
  plot = get_plot("heatmap", wildcards.plot_id)
  scaling = plot.get("harmonize_scaling", "NONE")
  if scaling not in ("CONTRASTS", "OVERALL"):
    return []
  return [f"results/jacusa2/cond1~{c['cond1']}/cond2~{c['cond2']}/bam~{wildcards.bam_type}/permutation.tsv"
          for c in pep.config["qutrna2"]["contrasts"]
          if scaling == "OVERALL" or plot in heatmap_plots(c["cond1"], c["cond2"])]


def _normalized_set_input(wildcards):
  return {**_plot_normalized_heatmap_input(wildcards), "trna_annotation": TRNA_ANNOTATION,
          "scale_with": _scale_tables(wildcards)}


def norm_score_label(plot):
  """The error types summed in an entry's normalized heatmap, e.g. mismatch + deletion."""
  columns = plot.get("norm_score", DEFAULT_NORM_SCORE).split("+")
  return " + ".join(c.removeprefix("perm_").removesuffix("_score_minus_max_null") for c in columns)


# each plots.heatmap entry again, with the normalized score of the permutation test
rule plot_normalized_heatmap_set:
  input: unpack(_normalized_set_input)
  output: pdf="results/plots/normalized/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}/normalized_heatmap.pdf",
          png="results/plots/normalized/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}/normalized_heatmap.png"
  wildcard_constraints: plot_id="[^/~]+"
  conda: "qutrna2"
  log: "logs/plot/normalized_heatmap/cond1~{COND1}/cond2~{COND2}/{plot_id}/bam~{bam_type}.log"
  params: sprinzl=lambda wildcards, input: f"--sprinzl {input.sprinzl}" if "sprinzl" in input.keys() else "",
          fdr=config["jacusa2"]["fdr"],
          opts=_normalized_set_opts,
          linker5=pep.config["qutrna2"]["linker5"],
          title=lambda wildcards: f"{wildcards.COND1}, {wildcards.COND2}, normalized score, "
                                  f"{norm_score_label(get_plot('heatmap', wildcards.plot_id))}"
  shell: """
    python {workflow.basedir:q}/scripts/plot_norm_heatmap.py \
        --conditions '{wildcards.COND1},{wildcards.COND2}' \
        --fdr {params.fdr} \
        --five_adapter {params.linker5} \
        --title {params.title:q} \
        {params.opts} \
        --out {output.pdf:q} \
        {params.sprinzl} \
        {input.perm:q} > {log:q} 2>&1
  """


rule plot_record_count:
  input: "results/stats/record_count.txt"
  output: pdf="results/plots/record_count/{FEATURE}.pdf",
          png="results/plots/record_count/{FEATURE}.png",
          rds="results/plots/record_count/{FEATURE}.rds",
  conda: "qutrna2"
  log: "logs/plot/record_count/{FEATURE}.log"
  params: basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         --type {wildcards.FEATURE} \
         --output {output.pdf:q} {input:q} \
         2> {log:q}
  """
rule plot_record_count_custom:
  input: "results/stats/record_count.txt"
  output: pdf="results/plots/record_count/custom/{plot_id}.pdf",
          png="results/plots/record_count/custom/{plot_id}.png",
          rds="results/plots/record_count/custom/{plot_id}.rds"
  conda: "qutrna2"
  log: "logs/plot/read_record_custom/{plot_id}.log"
  params: opts= lambda wildcards: get_plot("record_count", wildcards.plot_id).get("opts",""),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_record_count.R \
         {params.opts} \
         --output {output.pdf:q} {input:q} \
         2> {log:q}
  """


rule plot_read_length:
  input: "results/stats/read_length.txt"
  output: pdf="results/plots/read_length/{FEATURE}.pdf",
          png="results/plots/read_length/{FEATURE}.png",
          rds="results/plots/read_length/{FEATURE}.rds"
  conda: "qutrna2"
  log: "logs/plot/read_length/{FEATURE}.log"
  params: basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         --type {wildcards.FEATURE} \
         --output {output.pdf:q} {input:q} \
         2> {log:q}
  """
rule plot_read_length_custom:
  input: "results/stats/read_length.txt"
  output: pdf="results/plots/read_length/custom/{plot_id}.pdf",
          png="results/plots/read_length/custom/{plot_id}.png",
          rds="results/plots/read_length/custom/{plot_id}.rds"
  conda: "qutrna2"
  log: "logs/plot/read_length_custom/{plot_id}.log"
  params: opts=lambda wildcards: get_plot("read_length", wildcards.plot_id).get("opts", ""),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_read_length.R \
         {params.opts} \
         --output {output.pdf:q} {input:q} \
         2> {log:q}
  """


rule plot_threshold_summary:
  input: score="results/stats/alignment_score.txt",
         cutoff="results/stats/cutoff.txt"
  output: pdf="results/plots/alignment/threshold_summary.pdf",
          png="results/plots/alignment/threshold_summary.png",
          rds="results/plots/alignment/threshold_summary.rds"
  conda: "qutrna2"
  log: "logs/plot/threshold_summary.log"
  params: bam_types=",".join(["mapped", "mapped-random"]),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_threshold_summary.R \
         --type {params.bam_types:q} \
         --cutoff {input.cutoff:q} \
         --output {output.pdf:q} \
         {input.score:q} \
         2> {log:q}
  """
rule plot_threshold_summary_custom:
  input: score="results/stats/alignment_score.txt",
    cutoff="results/stats/cutoff.txt"
  output: pdf="results/plots/alignment/threshold_summary/custom/{plot_id}.pdf",
          png="results/plots/alignment/threshold_summary/custom/{plot_id}.png",
          rds="results/plots/alignment/threshold_summary/custom/{plot_id}.rds"
  conda: "qutrna2"
  log: "logs/plot/threshold_summary_custom/{plot_id}.log"
  params: opts=lambda wildcards: get_plot("threshold_summary", wildcards.plot_id).get("opts", ""),
          bam_types=",".join(["mapped", "mapped-random"]),
          basedir=workflow.basedir
  shell: """
    export QUTRNA2="{params.basedir}"
    Rscript {workflow.basedir}/scripts/plot_threshold_summary.R \
         {params.opts} \
         --type {params.bam_types:q} \
         --cutoff {input.cutoff:q} \
         --output {output.pdf:q} \
         {input.score:q} \
         2> {log:q}
  """
