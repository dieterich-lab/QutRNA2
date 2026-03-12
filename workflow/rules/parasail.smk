global REF_FASTA


########################################################################################################################
# Use GPU-assisted or only parasail to map reads

if config["alignment"]["method"] == "gpu":
  include: "parasail_gpu.smk"
else:
  if config["parasail"]["lines"] > 0:
    include: "parasail_split.smk"
  else:
    include: "parasail_map.smk"

########################################################################################################################
# filter alignments by random score distribution

rule parasail_infer_cutoff:
  input: real="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~real/{BC}_stats/alignment_score.txt",
         random="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/alignment~random/{BC}_stats/alignment_score.txt"
  output: score_plot="results/plots/alignment/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}/alignment_score.pdf",
          cutoff="results/bam/mapped/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}_stats/cutoff.tsv"
  # use one global precision target, then compute one cutoff per exact trna name.
  # this keeps the cutoff table directly aligned with observed references.
  params: precision=config["alignment"]["precision"]
  log: "logs/parasail/infer_cutoff/sample~{SAMPLE}/subsample~{SUBSAMPLE}/{BC}.log"
  conda: "qutrna2"

  shell: """
    python {workflow.basedir}/scripts/aln_score_cutoff_trna.py \
    --real {input.real:q} --random {input.random:q} \
    --precision {params.precision:q} \
    --score-plot {output.score_plot:q} \
    --output {output.cutoff:q} 2> {log:q}
  """
