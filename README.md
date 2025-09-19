# QutRNA2

Robust tRNA modification discovery from Nanopore direct tRNA sequencing

## New Features

QutRNA2 features the novel GPU-assisted [gpu-tRNA-mapper](https://github.com/fkallen/gpu-tRNA-mapper) that performs up to 25x faster than the previously used mapper [parasail](https://github.com/jeffdaily/parasail) for the same task. Furthermore, a new, improved version of [JACUSA v2.1.15](https://github.com/dieterich-lab/JACUSA2/releases/download/v2.1.15-RC/JACUSA_v2.1.15-RC.jar) is included that features subsampled scores that improve the signal-to-noise ratio when identifying tRNA modifications.

Finally, a filter framework has been added to the analysis workflow to remove spurious alignments by applying the following filters:

* Filter Random alignments
* Adapter overlap (with 5' and 3' splint adapters)
* Filter multimapper

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
We recommend the following installation order:

1. [QutRNA](https://github.com/dieterich-lab/QutRNA2)
2. [JACUSA2 2.1.15](https://github.com/dieterich-lab/JACUSA2/releases/download/v2.1.15-RC/JACUSA_v2.1.15-RC.jar)
3. [parasail v2.6.2](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz)
4. [gpu-tRNA-mapper](https://github.com/fkallen/gpu-tRNA-mapper)

### QutRNA2
Clone the repository and install the requirements with [conda](https://docs.conda.io/en/latest/).

Go to your desired `<QUTRNA2-LOCAL-DIR>` and clone the repository:
```console
cd <QUTRNA2-LOCAL-DIR>
git clone https://github.com/dieterich-lab/QutRNA2
```

Next, install all the requirements: 
```console
cd QutRNA2
conda env create -f conda-lock.yaml -n qutrna2
```
Finally, activate the environment:
```console
conda activate qutrna2
```
Make sure to activate the environment when you build the other dependencies!

### JACUSA2
Download [JACUSA2 2.1.15]([https://github.com/dieterich-lab/JACUSA2TODO](https://github.com/dieterich-lab/JACUSA2/releases/download/v2.1.15-RC/JACUSA_v2.1.15-RC.jar)) and store the JAR on your local computer. Put the JAR into a directory in your PATH or remember the location <JACUSA2_LOCAL_JAR> for later configuration.

### parasail
First, download the source code from parasail [v2.6.2.tar.gz](https://github.com/jeffdaily/parasail/archive/refs/tags/v2.6.2.tar.gz), store it in <PARASAIL-LOCAL-DIR>, and unpack it.
```console
cd <PARASAIL-LOCAL-DIR>
tar xzvpf v2.6.2.tar.gz
```
The `<PARASAIL-LOCAL-DIR>/parasail-2.6.2` directory is essential for later use. It contains the source code of [parasail](https://github.com/jeffdaily/parasail/).

Now, you can build [parasail](https://github.com/jeffdaily/parasail/):
```console
cd <PARASAIL-LOCAL-DIR>/parasail-2.6.2
mkdir build
cd build
cmake .. 
make
cd ../..
```

### gpu-tRNA-mapper
Install the [Software Requirements](https://github.com/fkallen/gpu-tRNA-mapper?tab=readme-ov-file#software-requirements) before building the [gpu-tRNA-mapper](https://github.com/fkallen/gpu-tRNA-mapper).

Once you have installed a CUDA toolkit and C++17 compiler, you can clone the [gpu-tRNA-mapper](https://github.com/fkallen/gpu-tRNA-mapper) repository:
```console
cd <GPU-TRNA-MAPPER-LOCAL-DIR>
git clone --recurse-submodules git@github.com:fkallen/gpu-tRNA-mapper.git
```

Check the [Setup](https://github.com/fkallen/gpu-tRNA-mapper?tab=readme-ov-file#setup) for details on configuring and building from source.

When you are in `<GPU-TRNA-MAPPER-LOCAL-DIR>`, you can start building the gpu-tRNA-mapper - use `<PARASAIL-LOCAL-DIR>` from building parasail:
```console
make \
	PARASAIL_INCLUDE_DIR=<PARASAIL-LOCAL-DIR>/parasail-2.6.2/include \
	PARASAIL_LIB_DIR=<PARASAIL-LOCAL-DIR>/parasail-2.6.2/lib \
	TARGET_GPU_ARCH=all 
```

Finally, install the mapper to `<GPU-TRNA-MAPPER-PATH>`:
```console
make install <GPU-TRNA-MAPPER-PATH>/bin
```

The binary is called `gpu-tRNA-mapper`.


## Setup QutRNA2 analysis
QutRNA2 uses YAML files to define the data (data.yaml) and parametrize the analysis (analysis.yaml). Finally, a TSV file provides the sample description.

In summary, the sample description `<SAMPLE_DESC>` must be a TAB-separated file and contain the following columns:

| condition | sample_name | subsample_name | base_calling | fastq\|bam |
| --------- | ----------- | -------------- | ------------ | ---------- |
| ...       | ...         | ...            | ...          | ...        |

See the files in the `QutRNA2/examples` folder. There is documentation for the available options.

## Setup data configuration
First, define your `<SAMPLE_DESC>`. This file holds sample-specific information, such as "condition", "sample_name", "subsample", and "fastq" - they directly correspond to columns - see `examples/sample_desc.tsv`.
Data for entries with the same "sample_name" will be merged - they represent technical replicates. For historical reasons, the column "base_calling" is present. Set it to "pass". Finally, the column "fastq" should point to the path of the reads corresponding to the specified entry.

Second, define your `<DATA_YAML>`. This file describes what reference and Sprinzl coordinates (if any) to use. See `examples/data.yaml`. Make sure to add your `<SAMPLE_DESC>`. Provide "ref_fasta" and define what Sprinzl coordinates to use and the size of the adapters used! Correct adapter lengths are essential!

### Sprinzl
For eukaryotic nuclear tRNAs, we use the following covariance model [TRNAinf-euk.cm](https://github.com/UCSC-LoweLab/tRAX/blob/master/TRNAinf-euk.cm) and labeling `data/nuclear-euk-masked.txt`.

For human mt-tRNAs, we use the sequence to Sprinzl mapping in [https://www.nature.com/articles/s41467-020-18068-6](https://www.nature.com/articles/s41467-020-18068-6) 
and deposited the data along with the Sprinzl labels to: `data/human_mt_seq_to_sprinzl.tsv` and `data/human_mt_sprinzl_labels.txt`. If you provide labels, ensure the first label is a "-"!

It is crucial to obtain covariance models for the organism and tRNAs studied. These models can be acquired, for example, from [https://github.com/UCSC-LoweLab/tRNAscan-SE/tree/master/lib/models](https://github.com/UCSC-LoweLab/tRNAscan-SE/tree/master/lib/models).

## Setup analysis configuration
Finally, define `<ANALYSIS_YAML>`. Here, the workflow is manipulated, and custom plots are defined. Check `examples/analysis.yaml` for examples. Add any necessary init code for GPU and provide paths for JACUSA2 and the gpu-tRNA-mapper, if they are unavailable in the standard path.

## Execute workflow
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
        --directory <ANALYSIS_OUTPUT>
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

## Results
When the analysis is finished, the `<ANALYSIS_OUTPUT>` directory will contain the following subdirectories:
"data", "info", "logs", and "results".

`<ANALYSIS_OUTPUT>/data/` will contain all unprocessed data used in the analysis.
`<ANALYSIS_OUTPUT>/info/` will contain some runtime information and the used configuration files to track used parameters.
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


