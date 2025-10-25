import os

import pandas as pd
from snakemake.io import temp, unpack
import re
import subprocess
import yaml


rule report_render:
  input: report="results/report/report.Rmd",
         preprocessing="results/report/preprocessing.Rmd",
         results="results/report/results.Rmd",
         data_config="results/report/data_config.Rmd",
         params="results/report/params.yaml"
  output: "results/report/report.html"
  log: "logs/report/render.log"
  shell: """
    mkdir -p results/report
    Rscript {workflow.basedir}/scripts/report_render.R \
      --params {input.params:q} \
      --format html_document \
      {input.report:q} \
      > {output:q} \
      2> {log:q}
  """


def _report_parse_template(wildcards):
  d = {
    "params": "results/report/params.yaml",
    "template": f"{workflow.basedir}/report/{{template}}.jinja"
  }
  if wildcards.template != "report":
    d["template"] = temp(d["template"])

  # TODO add dependency

  return d


rule report_parse_template:
  input: unpack(_report_parse_template)
  output: "results/report/{template}.Rmd"
  log: "logs/report_parse_template/{template}.log"
  conda: "qutrna2"
  shell: """
    mkdir -p results/report
    python {workflow.basedir}/scripts/report_parse_template.py \
      --output {output:q} \
      --params {input.params:q} \
      {input.template:q} \
      2> {log:q}
  """


# gather plots


def _report_create_params(wildcards):
  targets = [config["pepfile"], ]
  for fname in workflow.configfiles:
    targets.append(str(fname))

  # add feature rds
  targets.extend(flatten_dict(get_read_feature_plots("rds")))

  # add heatmap plots files.txt
  conds1 = [contrast["cond1"] for contrast in pep.config["qutrna2"]["contrasts"]]
  conds2 = [contrast["cond2"] for contrast in pep.config["qutrna2"]["contrasts"]]
  plot_ids = [plot["id"] for plot  in config["plots"]["heatmap"]]
  bam_types = []
  if config["call_filtered"]:
    bam_types = [f"filtered-{f}" for f in FILTERS_APPLIED]
  bam_types.append("final")
  targets.extend(expand("results/plots/scores/cond1~{cond1}/cond2~{cond2}/{plot_id}/bam~{bam_type}/files.tsv",
    cond1=conds1, cond2=conds2,
    plot_id=plot_ids,
    bam_type=bam_types))

  return targets

rule report_create_params:
  input: _report_create_params
  output: "results/report/params.yaml"
  run:
    d = {
      "config": config,
      "config_str": yaml.dump(config),
      "pep": pep.to_dict(),
      "pep_str": yaml.dump(pep.to_dict()),
      "basedir": workflow.basedir,
      "workdir": os.getcwd(),
      "configfiles": [str(f) for f in workflow.configfiles],
      "pepfile": config["pepfile"],
      "filters_applied": FILTERS_APPLIED,
    }
    versions = {
      "qutrna2": VERSION,
    }
    d["rds"] = get_read_feature_plots("rds")
    d["rds"]["heatmap"] = get_heatmap_plots("rds")

    sample_table = pd.read_csv(pep.config["sample_table"],sep="\t")
    d["sample_table"] = {}
    d["sample_table"]["header"] = sample_table.columns.to_list()
    d["sample_table"]["rows"] = sample_table.to_records().tolist()

    try:
      jacusa2_return = subprocess.check_output([config["jacusa2"]["jar"], "-h"])
      m = re.search(r"Version:\t([^\n]+)", jacusa2_return.decode("utf-8"), re.M)
      versions["jacusa2"] = m.group(1)
    except Exception:
      versions["jacusa2"] = "unknown"

    try:
      infernal_return = subprocess.check_output(["cmalign", "-h"])
      m = re.search(r"\n# INFERNAL ([^\n]+)", infernal_return.decode("utf-8"), re.M)
      versions["infernal"] = m.group(1)
    except Exception:
      versions["infernal"] = "unknown"

    # FIXME gpu-tRNA-mapper does not support version information
    #try:
    #  gpu_trna_mapper_return = subprocess.check_output([config["gpu"]["bin"], "-h"])
    #  m = re.search(r"\n# INFERNAL ([^\n]+)", gpu_trna_mapper_return.decode("utf-8"), re.M)
    #  versions["gpu_trna_mapper"] = m.group(1)
    #except Exception:
    versions["gpu_trna_mapper"] = "unknown"

    d["versions"] = versions

    try:
      os.mkdir("results/report")
    except FileExistsError:
      pass
    with open(output[0], "w") as f:
      yaml.safe_dump(d, f)
