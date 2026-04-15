# QutRNA2

Robust tRNA modification discovery from Nanopore direct tRNA sequencing

If you use QutRNA2, please cite: [https://www.biorxiv.org/content/10.1101/2025.10.20.683443v1](https://www.biorxiv.org/content/10.1101/2025.10.20.683443v1)

## Table of Contents

- [Quick Start](#quick-start)
- [New Features](#new-features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Conda](#environment-setup-using-conda)
  - [Singularity / Apptainer](#using-singularity-container)
- [Setup QutRNA2 Analysis](#setup-qutrna2-analysis)
- [Setup Data Configuration](#setup-data-configuration)
  - [Sprinzl](#sprinzl)
- [Setup Analysis Configuration](#setup-analysis-configuration)
- [Examples](#examples)
- [Execute Workflow](#execute-workflow)
- [Results](#results)
  - [Alignments](#alignments)
  - [cmalign](#cmalign)
  - [JACUSA2](#jacusa2)
  - [Plots](#plots)
  - [Stats](#stats)
  - [Secondary Structure](#secondary-structure-ss)

## Quick Start

This quick start guide assumes you have [conda](https://docs.conda.io/en/latest/) installed and a CUDA-capable GPU available. See [Installation](#installation) for full details and a container solution using Singularity.

**1. Clone the repository:**
```console
git clone https://github.com/dieterich-lab/QutRNA2
cd QutRNA2
```

**2. Install and activate the environment:**
```console
conda env create -f conda.yaml -n qutrna2
conda activate qutrna2
```

**3. Copy and edit the example config files:**

```console
cp examples/analysis/map_with_gpu.yaml my_analysis.yaml
cp examples/data/sprinzl_cm.yaml my_data.yaml
```

Edit `my_data.yaml` to point to your reference FASTA [examples](https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/Hsapi19-seq.html), sample description TSV (labels and filepaths to basecalled samples) [example](https://github.com/dieterich-lab/QutRNA2/blob/main/examples/sample_table_fastq.tsv), and Sprinzl coordinate labels [example for eukaryotic cytosolic tRNA](https://github.com/dieterich-lab/QutRNA2/blob/main/data/nuclear-euk-masked.txt).
Edit `my_analysis.yaml` to set paths for JACUSA2 and gpu-tRNA-mapper, and any GPU init commands.

**4. Run a dry run to verify the workflow:**
```console
snakemake \
    -c 1 \
    --snakefile <QUTRNA2_LOCAL_DIR>/workflow/Snakefile \
    --use-conda \
    --configfiles my_analysis.yaml \
    --config pepfile=my_data.yaml \
    --directory <ANALYSIS_OUTPUT> \
    -n
```

A list of jobs to be executed should appear. If it does without errors, remove `-n` and increase `-c` to the number of available cores to start the analysis.

## New Features

QutRNA2 features the novel GPU-assisted [gpu-tRNA-mapper](https://github.com/fkallen/gpu-tRNA-mapper) that performs up to 25x faster than the previously used mapper [parasail](https://github.com/jeffdaily/parasail) for the same task. Furthermore, a new, improved version of [JACUSA v2.1.16](https://github.com/dieterich-lab/JACUSA2/releases/download/v2.1.16/JACUSA_v2.1.16.jar) is included, featuring subsampled scores that improve the signal-to-noise ratio for identifying tRNA modifications.

Finally, a filter framework has been added to the analysis workflow to remove spurious alignments by applying the following filters:

* Filter Random alignments
* Filter Adapter overlap (with 5' and 3' splint adapters)
* Filter multimapping reads

We added the following plots to assess the impact of filtering:

* Alignment threshold summary plot
* Impact of filters on read length
* Impact of filters on the number of reads

More customization options for heatmap plots:

* Filter tRNAs by min. number of reads
* Display or ignore specific tRNAs by regular expression
* Mark positions of interest
* Use patterns to customize the title of heatmap plots

## Requirements
To use GPU-assisted mapping, you need a compatible NVIDIA GPU. For details, check [Hardware requirements](https://github.com/fkallen/gpu-tRNA-mapper?tab=readme-ov-file#hardware-requirements).
In brief, a CUDA-capable GPU with Volta architecture or newer is recommended.

If no compatible GPU is present, QutRNA2 can be used with [parasail](https://github.com/jeffdaily/parasail) but will run much slower.

## Installation

### Environment setup using Conda

We provide a conda file with all necessary packages. 
Clone the repository and install the requirements with [conda](https://docs.conda.io/en/latest/).

Go to your desired `<QUTRNA2-LOCAL-DIR>` and clone the repository:
```console
cd <QUTRNA2-LOCAL-DIR>
git clone https://github.com/dieterich-lab/QutRNA2
```

Next, install all the requirements: 
```console
cd QutRNA2
conda env create -f conda.yaml -n qutrna2
```

Finally, activate the environment:
```console
conda activate qutrna2
```
### Using Singularity container

If Singularity/Apptainer is not already available on your system, see the [SingularityCE](https://docs.sylabs.io/guides/latest/admin-guide/installation.html) or [Apptainer](https://apptainer.org/docs/admin/latest/installation.html) installation guides.

A Singularity definition file is provided at [`singularity/qutrna2.def`](https://github.com/dieterich-lab/QutRNA2/blob/main/singularity/qutrna2.def) to build a portable container image.

**Build the image:**
```console
singularity build qutrna2.sif singularity/qutrna2.def
```

> Note: building a Singularity image typically requires root or `--fakeroot` privileges (e.g. `singularity build --fakeroot qutrna2.sif singularity/qutrna2.def`)
> If building on an HPC cluster, check whether `--fakeroot` is supported.

## Setup QutRNA2 analysis
QutRNA2 uses YAML files to define the data (data.yaml) and parametrize the analysis (analysis.yaml). Finally, a TSV file provides the sample description.

In summary, the sample description `<SAMPLE_DESC>` must be a TAB-separated file and contain the following columns:

| condition | sample_name | subsample_name | base_calling | fastq\|bam |
| --------- | ----------- | -------------- | ------------ | ---------- |
| ...       | ...         | ...            | ...          | ...        |

See the files in the `QutRNA2/examples` folder for documented YAML and toy examples for sample tables. QutRNA2 distinguishes the configuration of the analysis and the data. The following analysis types are supported: 

* map reads with gpu-tRNA-mapper (see `QutRNA2/examples/analysis/map_with_gpu.yaml`),
* map reads with parasail (see `QutRNA2/examples/analysis/map_with_parasail.yaml`), and
* use exisiting mapping (see `QutRNA2/examples/analysis/existing_mapping.yaml`.

QutRNA2 supports the following approaches to assign Sprinzl coordinates :

* using a covarince model and secondary structure alignment (see `QutRNA2/examples/data/sprinzl_cm.yaml`),
* using an existing aligned FASTA file (see `QutRNA2/examples/data/sprinzl_afasta.yaml`), or
* a direct mapping of sequence to Sprinzl coordinates (see `QutRNA2/examples/data/seq_to_sprinzl.yaml`)).

Those files are templates and must be adjusted to the user's needs.

## Setup data configuration
First, define your `<SAMPLE_DESC>`. This file holds sample-specific information, such as "condition", "sample_name", "subsample", and "fastq" or "bam" - they directly correspond to columns - see `examples/sample_desc_fastq.tsv`.
Data for entries with the same "sample_name" will be merged - they represent technical replicates. For historical reasons, the column "base_calling" is present. Set it to "pass". Finally, the column "fastq" should point to the path of the compressed (gzip) fastq file.

Second, define your `<DATA_YAML>`. This file describes what reference and Sprinzl coordinates (if any) to use. See `examples/data/*.yaml`. Make sure to add your `<SAMPLE_DESC>`. Provide "ref_fasta" and define what Sprinzl coordinates to use and the size of the adapters used! Correct adapter lengths are essential!

### Sprinzl
For eukaryotic nuclear tRNAs, we use the following covariance model [TRNAinf-euk.cm](https://github.com/UCSC-LoweLab/tRAX/blob/master/TRNAinf-euk.cm) and labeling `data/nuclear-euk-masked.txt`.

For human mt-tRNAs, we use the sequence to Sprinzl mapping in [https://www.nature.com/articles/s41467-020-18068-6](https://www.nature.com/articles/s41467-020-18068-6) 
and deposited the data along with the Sprinzl labels to: `data/human_mt_seq_to_sprinzl.tsv` and `data/human_mt_sprinzl_labels.txt`.

It is crucial to obtain covariance models for the organism and tRNAs studied. These models can be acquired, for example, from [https://github.com/UCSC-LoweLab/tRNAscan-SE/tree/master/lib/models](https://github.com/UCSC-LoweLab/tRNAscan-SE/tree/master/lib/models).

## Setup analysis configuration
Finally, define `<ANALYSIS_YAML>`. Here, the workflow is manipulated, and custom plots are defined. Check `examples/analysis/*.yaml` for examples. Add any necessary init code for the GPU and provide paths for JACUSA2 and gpu-tRNA-mapper if they are not in the standard path.

## Examples

## Execute workflow

### Run using the Conda environment
If not done yet, activate qutrna2 conda environment:
```console
conda activate qutrna2
```

Use `<ANALYSIS_OUTPUT>` folder to define where QutRNA2 should write the output to:
```console
snakemake \
    -c 1 \
    --snakefile <QUTRNA_LOCAL_DIR>/Snakefile \
    --use-conda \
    --configfiles <ANALYSIS_YAML> \
        --config pepfile=<DATA_YAML> \
        --directory <ANALYSIS_OUTPUT> \
    -n 
```
You should see a list of necessary jobs to be run.
You should increase "-c 1" to whatever suits your computing machine.

Now, you can start the analysis (remove "-n")
```console
snakemake \
    -c 1 \
    --snakefile <QUTRNA_LOCAL_DIR>/Snakefile \
    --use-conda \
    --configfiles <ANALYSIS_YAML> \
        --config pepfile=<DATA_YAML> \
        --directory <ANALYSIS_OUTPUT>
```

### Run using the container

The Snakemake command is now the same as above, except it is called as part of a Singularity command:

```console
singularity run \
  --nv \
  --bind <HOST_DATA_DIR>:<HOST_DATA_DIR> \
  <PATH_TO_SIF>/qutrna2.sif \
  snakemake \
    --cores <N_CORES> \
    --snakefile <QUTRNA2_LOCAL_DIR>/workflow/Snakefile \
    --configfiles <ANALYSIS_YAML> \
        --config pepfile=<DATA_YAML> \
        --directory <ANALYSIS_OUTPUT>
```

- `singularity run --nv` enables GPU access inside the container (required for gpu-tRNA-mapper).
- `singularity run --bind` mounts a host directory into the container so input/output paths remain accessible. Bind all directories referenced in your YAML configs or a parent directory.
- `snakemake --executor local` is recommended when running inside a SLURM job (i.e. resource allocation is already handled by SLURM). Set `--cores` to the number of cores allocated to your job (e.g. `$SLURM_CPUS_ON_NODE`).

## Results
When the analysis is finished, the `<ANALYSIS_OUTPUT>` directory will contain the following subdirectories:
"data", "info", "logs", and "results".

`<ANALYSIS_OUTPUT>/data/` will contain all unprocessed data used in the analysis.
`<ANALYSIS_OUTPUT>/info/` will contain runtime information and the configuration files used to track parameters.
`<ANALYSIS_OUTPUT>/logs/` will contain logs for executed jobs.
`<ANALYSIS_OUTPUT>/results/` will contain all calculations.

### data
The directory `<ANALYSIS_OUTPUT>/results/data` will contain processed instances of the reference sequence.

### Alignments
Alignments are stored in `<ANALYSIS_OUTPUT>/results/bam/<read-type>/...`. `<read-type>` corresponds to mapped, filtered, and final reads. The BAMs in the subdirectory "final" are used to calculate JACUSA2 score profiles.
Each subdirectory is organised according to "sample_name", "subsample_name", and "base_calling" columns from `<SAMPLE_DESC>`.

### cmalign
If a covariance model was provided in `<DATA_YAML>`, the secondary structure alignment under `<ANALYSIS_OUTPUT>/results/cmalign/align.stk` will be generated.

### jacusa2
The directory `<ANALYSIS_OUTPUT>/results/jacusa2` will contain JACUSA2 results for defined contrasts.

## plots
The directory `<ANALYSIS_OUTPUT>/results/plots` will contain plots.

Check `<ANALYSIS_OUTPUT>/results/plots/scores/cond1~{cond1}/cond2~{cond1}/{id}/bam~final/heatmap.pdf` for JACUAS2 score profiles after filtering. This plot concludes the analysis.

### stats
If filters were applied, the directory `<ANALYSIS_OUTPUT>/results/stats` will contain summary statistics for features such as alignment score, read length, and read count.

### secondary structure (ss)
`<ANALYSIS_OUTPUT>/results/seq_to_sprinzl_filtered.tsv` will contain the sequence to sprinzl mapping.

