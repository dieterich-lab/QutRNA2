import pandas as pd

global create_include
global REF_FILTERED_TRNAS_FASTA


########################################################################################################################
# Coordinate system: seq(uence) or sprinzl

# final sequence to sprinzl mapping
SEQ_TO_SPRINZL_FINAL = "results/ss/seq_to_sprinzl_filtered.tsv"

# in case of sprinzl, handle if model or precalculated mapping is to be used
if pep.config["qutrna2"]["coords"] == "sprinzl":
  # How mapping to sprinzl coordinates are calculated, by alignment or from existing mapping
  # cm -> covarince model -> calculate secondary structure alignment and create mapping
  # seq_to_sprinzl -> use exisiting mapping
  if "cm" in pep.config["qutrna2"]["sprinzl"] or "scheme" in pep.config["qutrna2"]["sprinzl"]:
    _SPRINZL_CFG = pep.config["qutrna2"]["sprinzl"]
    _VALID_SCHEMES = {"euk", "arch", "bact", "mito"}
    _HAS_CM = "cm" in _SPRINZL_CFG
    _HAS_SPRINZL_SCHEME = "scheme" in _SPRINZL_CFG

    if _HAS_CM and _HAS_SPRINZL_SCHEME:
      raise ValueError(
          "In CM mode, provide only one of qutrna2.sprinzl.scheme or qutrna2.sprinzl.cm."
      )
    if not _HAS_CM and not _HAS_SPRINZL_SCHEME:
      raise ValueError(
          "In CM mode, you must provide one of qutrna2.sprinzl.scheme or qutrna2.sprinzl.cm."
      )

    SPRINZL_MODE = "cm"
    SEQ_TO_SPRINZL_INIT = "results/ss/seq_to_sprinzl.tsv"

    # qutrna2.sprinzl.labels sets the ordered list of Sprinzl labels shown in
    # the heatmap; blank a position with '-' to drop it from the plot.
    # Required with cm (the only source of a label list); optional with
    # scheme, where sprinx assigns labels and this only curates the plot.
    if "labels" in _SPRINZL_CFG:
      SPRINZL_LABELS = "data/sprinzl_labels.txt"
      create_include("sprinzl_labels",
          _SPRINZL_CFG["labels"],
          SPRINZL_LABELS,
          config["include"].get("sprinzl_labels", "copy"))
    elif _HAS_CM:
      raise ValueError("qutrna2.sprinzl.cm requires qutrna2.sprinzl.labels.")

    if _HAS_SPRINZL_SCHEME:
      # sprinx runs its own cmalign internally and assigns Sprinzl labels
      # from stem-loop topology, so this branch skips the AFASTA/consensus-
      # label chain entirely.
      _SPRINZL_SCHEME = _SPRINZL_CFG["scheme"]
      if _SPRINZL_SCHEME not in _VALID_SCHEMES:
        raise ValueError(
            f"Invalid qutrna2.sprinzl.scheme='{_SPRINZL_SCHEME}'. "
            f"Must be one of: {', '.join(sorted(_VALID_SCHEMES))}"
        )

      _MITO_CFG = _SPRINZL_CFG.get("mito", {})
      _CYTO_CFG = _SPRINZL_CFG.get("cyto", {})
      if _SPRINZL_SCHEME == "mito" and _CYTO_CFG:
        raise ValueError("qutrna2.sprinzl.cyto only applies when scheme is euk/arch/bact.")
      if _SPRINZL_SCHEME != "mito" and _MITO_CFG:
        raise ValueError("qutrna2.sprinzl.mito only applies when scheme=mito.")

      _SPRINX_OUT = "results/ss/sprinx_mapping.tsv"
      _sprinx_opts = ["--debug"]
      if "canonical_cm" in _MITO_CFG:
        _sprinx_opts.append("--canonical-cm " + " ".join(_MITO_CFG["canonical_cm"]))
      if "armless_cm_dir" in _MITO_CFG:
        _sprinx_opts.append(f"--armless-cm-dir {_MITO_CFG['armless_cm_dir']}")
      if "cyto_cm_db" in _CYTO_CFG:
        _sprinx_opts.append(f"--cyto-cm-db {_CYTO_CFG['cyto_cm_db']}")
      _SPRINX_OPTS = " ".join(_sprinx_opts)

      rule ss_check_sprinx_headers:
        input: REF_FILTERED_TRNAS_FASTA
        output: touch("results/ss/sprinx_headers_ok")
        log: "logs/ss/check_sprinx_headers.log"
        conda: "qutrna2"
        shell: """
          python {workflow.basedir}/scripts/sprinzl_utils.py check-sprinx-headers \
            {input:q} \
            2> {log:q}
        """

      rule ss_sprinx_run:
        input: fasta=REF_FILTERED_TRNAS_FASTA,
               headers_ok="results/ss/sprinx_headers_ok"
        output: _SPRINX_OUT
        log: "logs/ss/sprinx_run.log"
        params: scheme=_SPRINZL_SCHEME, opts=_SPRINX_OPTS
        threads: workflow.cores
        conda: "qutrna2"
        shell: """
          sprinx --scheme {params.scheme} \
            --fasta {input.fasta:q} \
            --out {output:q} \
            -p {threads} \
            {params.opts} \
            2> {log:q}
        """

      rule ss_sprinx_to_seq_to_sprinzl:
        input: _SPRINX_OUT
        output: SEQ_TO_SPRINZL_INIT
        conda: "qutrna2"
        log: "logs/ss/sprinx_to_seq_to_sprinzl.log"
        shell: """
          python {workflow.basedir}/scripts/sprinzl_utils.py sprinx-to-seq-to-sprinzl \
            --output {output:q} \
            {input:q} \
            2> {log:q}
        """
    else:
      # covariance model destination
      CM = "data/cm.stk"
      AFASTA = "results/ss/ref.afasta"
      create_include("cm",
        pep.config["qutrna2"]["sprinzl"]["cm"],
        CM,
        config["include"]["cm"])
      CONSENSUS_LABELS = "results/data/consensus_labels.tsv"

      rule cmalign_run:
        input: cm=CM,
          fasta=REF_FILTERED_TRNAS_FASTA
        output: "results/cmalign/align.stk"
        log: "logs/cmalign/run.log"
        threads: 1
        params: opts=config["cmalign"]["opts"]
        conda: "qutrna2"
        shell: """
          cmalign {params.opts} \
            --cpu {threads} \
            -o {output:q} \
            {input.cm:q} \
            {input.fasta:q} \
            2> {log}
       """

      rule ss_stk_to_afasta:
        input: stk="results/cmalign/align.stk",
        output: AFASTA
        conda: "qutrna2"
        log: "logs/ss/stk_to_afasta.log"
        shell: """
          python {workflow.basedir}/scripts/sprinzl_utils.py stk-to-afasta \
            --output {output:q} \
            {input.stk:q} \
            2> {log:q}
        """

      rule ss_create_consensus_label:
        input: stk="results/cmalign/align.stk",
               labels=SPRINZL_LABELS,
        output: CONSENSUS_LABELS,
        log: "logs/ss/create_consus_labels.log"
        shell: """
          python {workflow.basedir}/scripts/sprinzl_utils.py consensus-labels \
          --labels={input.labels:q} \
          --output {output:q} \
          {input.stk:q} \
          2> {log:q}
        """

      rule ss_afasta_to_sprinzl:
        input: afasta=AFASTA,
               labels=CONSENSUS_LABELS
        output: SEQ_TO_SPRINZL_INIT
        conda: "qutrna2"
        log: "logs/ss/afasta_to_sprinzl.log"
        shell: """
          python {workflow.basedir}/scripts/sprinzl_utils.py afasta-to-sprinzl \
            --output {output:q} \
            --consensus-labels {input.labels:q} \
            {input.afasta:q} \
            2> {log:q}
        """
  elif "afasta" in pep.config["qutrna2"]["sprinzl"]:
    SPRINZL_MODE = "afasta"
    AFASTA = "data/ref.afasta"
    create_include("afasta",
      pep.config["qutrna2"]["sprinzl"]["afasta"],
      AFASTA,
      config["include"].get("afasta", "copy"))
    SEQ_TO_SPRINZL_INIT = "results/ss/seq_to_sprinzl.tsv"
    _CONSENSUS_LABELS = "data/consensus_labels.tsv"
    create_include("consensus_labels",
      pep.config["qutrna2"]["sprinzl"]["consensus_labels"],
      _CONSENSUS_LABELS,
      config["include"].get("consensus_labels", "copy"))
    SPRINZL_LABELS = _CONSENSUS_LABELS
    CONSENSUS_LABELS = "results/data/consensus_labels.tsv"

    rule process_consensus:
      input: _CONSENSUS_LABELS
      output: CONSENSUS_LABELS
      shell: """
        awk ' BEGIN {{ print "label" }} ; {{ print }} ' {input:q} > {output:q}
      """

    rule ss_afasta_to_sprinzl:
      input: afasta=AFASTA,
             labels=CONSENSUS_LABELS
      output: SEQ_TO_SPRINZL_INIT
      conda: "qutrna2"
      log: "logs/ss/afasta_to_sprinzl.log"
      shell: """
        python {workflow.basedir}/scripts/sprinzl_utils.py afasta-to-sprinzl \
          --output {output:q} \
          --consensus-labels {input.labels:q} \
          {input.afasta:q} \
          2> {log:q}
      """

  elif "seq_to_sprinzl" in pep.config["qutrna2"]["sprinzl"]:
    SPRINZL_MODE = "seq2sprinzl"
    SEQ_TO_SPRINZL_INIT = "data/seq_to_sprinzl.tsv"
    create_include("seq_to_sprinzl",
      pep.config["qutrna2"]["sprinzl"]["seq_to_sprinzl"],
      SEQ_TO_SPRINZL_INIT,
      config["include"].get("seq_to_sprinzl", "copy"))

    # optional override of the auto-derived plot label list; see the "cm"
    # branch above for what qutrna2.sprinzl.labels means.
    if "labels" in pep.config["qutrna2"]["sprinzl"]:
      SPRINZL_LABELS = "data/sprinzl_labels.txt"
      create_include("sprinzl_labels",
          pep.config["qutrna2"]["sprinzl"]["labels"],
          SPRINZL_LABELS,
          config["include"].get("sprinzl_labels", "copy"))
  else:
    raise Exception("Must provide either 'cm' or 'seq_to_sprinzl' in pep!")

  if "SPRINZL_LABELS" not in globals():
    # cm+scheme and seq_to_sprinzl both produce a seq_to_sprinzl TSV directly,
    # with no separately maintained ordered label list; derive the flat label
    # list plot_heatmap.R needs from the union of labels in it.
    SPRINZL_LABELS = "results/data/sprinzl_labels_auto.txt"

    rule ss_derive_sprinzl_labels:
      input: SEQ_TO_SPRINZL_FINAL
      output: SPRINZL_LABELS
      conda: "qutrna2"
      log: "logs/ss/derive_sprinzl_labels.log"
      shell: """
        python {workflow.basedir}/scripts/sprinzl_utils.py derive-sprinzl-labels \
          --output {output:q} \
          {input:q} \
          2> {log:q}
      """

else:
  SEQ_TO_SPRINZL_INIT = "data/seq_to_sprinzl.tsv"

  create_include("seq_to_sprinzl",
    pep.config["qutrna2"]["sprinzl"]["seq_to_sprinzl"],
    SEQ_TO_SPRINZL_INIT,
    config["include"].get("seq_to_sprinzl", "copy"))

########################################################################################################################
# Filter sprinzl mapping and apply

rule ss_seq_to_sprinzl_final:
  input: SEQ_TO_SPRINZL_INIT
  output: SEQ_TO_SPRINZL_FINAL
  run:
    df = pd.read_csv(input[0], sep="\t")
    df = df[df["sprinzl"] != "-"]
    df.to_csv(str(output[0]), sep="\t", index=False, quoting=False)


rule ss_transform:
  input: jacusa2="results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/scores_seq.tsv",
         seq_sprinzl=SEQ_TO_SPRINZL_FINAL
  output: "results/jacusa2/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}/scores_sprinzl.tsv"
  conda: "qutrna2"
  log: "logs/ss/transform/cond1~{COND1}/cond2~{COND2}/bam~{bam_type}.log"
  params: linker5=pep.config["qutrna2"]["linker5"]
  shell: """
    python {workflow.basedir:q}/scripts/sprinzl_utils.py transform \
        --sprinzl {input.seq_sprinzl:q} \
        --output {output:q} \
        --linker5 {params.linker5} \
        {input.jacusa2:q} 2> {log:q}
  """
