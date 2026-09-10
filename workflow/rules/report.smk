import glob
import importlib.metadata
import json
import os
import re
import shutil
import subprocess

import pandas as pd
import yaml


rule report_render:
  input: params="results/report/params.yaml",
         template=f"{workflow.basedir}/report/report.html.jinja"
  output: "results/report/report.html"
  log: "logs/report/render.log"
  conda: "qutrna2"
  shell: """
    python {workflow.basedir}/scripts/report_render.py \
      --params {input.params:q} \
      --output {output:q} \
      {input.template:q} \
      2> {log:q}
  """


# gather plots


def _report_create_params(wildcards):
  targets = [config["pepfile"], ]
  for fname in workflow.configfiles:
    targets.append(str(fname))

  # the report embeds these, so they have to exist before params are written
  targets.extend(flatten_dict(get_read_feature_plots("png")))

  # add heatmap plots files.txt
  bam_types = []
  if config["call_filtered"]:
    bam_types = [f"filtered-{f}" for f in FILTERS_APPLIED]
  bam_types.append("final")
  for contrast in pep.config["qutrna2"]["contrasts"]:
    cond1, cond2 = contrast["cond1"], contrast["cond2"]
    targets.extend(expand("results/plots/scores/cond1~{cond1}/cond2~{cond2}/{plot_id}/bam~{bam_type}/files.tsv",
      cond1=cond1, cond2=cond2,
      plot_id=[plot["id"] for plot in heatmap_plots(cond1, cond2)],
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
    # which checkout produced the run, so a report can be traced back to code
    if shutil.which("git"):
      def _git(*args):
        return subprocess.run(["git", *args], stdout=subprocess.PIPE,
                              cwd=workflow.basedir).stdout.decode().strip()
      rev = _git("rev-parse", "--short", "HEAD")
      if rev:
        versions["qutrna2"] += f" ({rev}{'+dirty' if _git('status', '--porcelain') else ''})"
    d["png"] = get_read_feature_plots("png")
    d["png"]["heatmap"] = get_heatmap_plots("png")

    bam_types = ["final"]
    if config["call_filtered"]:
      bam_types += [f"filtered-{f}" for f in FILTERS_APPLIED]

    # everything a reader might want to open themselves, listed whether or not
    # the report embeds it
    d["files"] = {
      "bam": sorted(expand("results/bam/final/{sample}.sorted.bam", sample=SAMPLES)),
      "scores": sorted(
        f"results/jacusa2/cond1~{c['cond1']}/cond2~{c['cond2']}/bam~{b}/{SCORES}"
        for c in pep.config["qutrna2"]["contrasts"] for b in bam_types),
      "plots": sorted(flatten_dict(d["png"])),
    }

    sample_table = pd.read_csv(pep.config["sample_table"],sep="\t")
    d["sample_table"] = {}
    d["sample_table"]["header"] = sample_table.columns.to_list()
    d["sample_table"]["rows"] = sample_table.to_records().tolist()

    def conda_version(package):
      """version recorded for an installed conda package, or None."""
      prefix = os.environ.get("CONDA_PREFIX")
      if not prefix:
        return None
      for fname in glob.glob(os.path.join(prefix, "conda-meta", f"{package}-*.json")):
        with open(fname) as f:
          return json.load(f).get("version")
      return None

    def pip_version(package):
      """installed version, plus the source commit when pip took it from a
      repository. A package pinned by revision can keep one version string
      across revisions, and then the commit is the only thing that identifies
      the code that ran."""
      try:
        version = importlib.metadata.version(package)
        url = importlib.metadata.distribution(package).read_text("direct_url.json")
        commit = json.loads(url or "{}").get("vcs_info", {}).get("commit_id", "")
        return f"{version} ({commit[:7]})" if commit else version
      except Exception:
        return None

    try:
      jacusa2_return = subprocess.check_output([config["jacusa2"]["jar"], "-h"])
      m = re.search(r"Version:\t([^\n]+)", jacusa2_return.decode("utf-8"), re.M)
      versions["jacusa2"] = m.group(1)
    except Exception:
      versions["jacusa2"] = conda_version("jacusa2") or "unknown"

    try:
      infernal_return = subprocess.check_output(["cmalign", "-h"])
      m = re.search(r"\n# INFERNAL ([^\n]+)", infernal_return.decode("utf-8"), re.M)
      versions["infernal"] = m.group(1)
    except Exception:
      versions["infernal"] = conda_version("infernal") or "unknown"

    # gpu-tRNA-mapper has no version flag, so the package metadata is the only
    # place a version is recorded
    versions["gpu_trna_mapper"] = conda_version("gpu-trna-mapper") or "unknown"
    versions["sprinx"] = pip_version("sprinx") or conda_version("sprinx") or "unknown"

    d["versions"] = versions

    try:
      os.mkdir("results/report")
    except FileExistsError:
      pass
    with open(output[0], "w") as f:
      yaml.safe_dump(d, f, sort_keys=False)
