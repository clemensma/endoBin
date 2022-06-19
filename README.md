> Note: endoMiner is in the construction stage 🚧.

# endoMiner

------------------

This pipeline separates endosymbiont prokaryotic genomes from their respective eukaryotic host mitogenome based on raw Illumina reads. It performs raw- and trimmed read quality control, adapter trimming, de-novo assembly, contig sorting and assembly quality assessment on paired short read sequence data.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.0-brightgreen.svg)](http://nextflow.io)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-MIT-lightgrey.svg)

## Installing nextflow

[Nextflow](https://www.nextflow.io/) is required to run this pipeline. For
installing Nextflow, simply follow the [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)
in their documentation.

## Choosing the right container service

This pipeline runs with [Docker](https://www.docker.com/) or [Singularity](https://singularity.hpcng.org/). If you want to run this pipeline locally with root privileges, we recommend using Docker. If this pipeline is executed on a high performance cluster, we recommend Singularity. Please ask your system administrator which is applicable.

## Downloading this pipeline

To download this pipeline, go to the folder where you wish to store this
pipeline and then run the following:

    $ git clone https://github.com/clemensma/endoMiner

## Quickstart

Launch the pipeline execution with the following command: 

    $ nextflow run main.nf -profile local

This will run the pipeline with default commands and default data. All required dependencies will be downloaded and installed in isolated containers.

## Order of execution:

1. Raw reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Trimmed reads quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. De Novo assembly ([SPAdes](https://cab.spbu.ru/software/spades/))
5. Contig sorting based on endosymbiont reference genome ([blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [bfg](https://github.com/fethalen/better_fasta_grep))
6. Read mapping for coverage estimate ([bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
7. Coverage estimate
8. Endosymbiont genome quality assessment ([BUSCO](https://busco.ezlab.org/), [CheckM](https://ecogenomics.github.io/CheckM/))
9. Mitogenome extraction ([blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [bfg](https://github.com/fethalen/better_fasta_grep))
10. Mitogenome reassembly ([NOVOPlasty](https://github.com/ndierckx/NOVOPlasty))
11. Mitogenome strand control ([blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi))
12. Mitogenome annotation ([MITOS](http://mitos.bioinf.uni-leipzig.de/index.py))
13. Formatting of MITOS output ([blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi))

## Input files

The endoMiner pipeline needs the following files:
* Paired-end Illumina reads, `*.fastq` or `*fastq.gz`
* A reference genome, `*.fa`

The Illumina read file names should match the following naming convention: *name_R{1,2}.extension*

where:
* *name* is the identifier of the sample;
* the number **1** or **2** is the read pair in the paired-end samples;
* *extension* is the read file name extension eg. `fq`, `fq.gz`, `fastq.gz`, etc.

example: `all_R{1,2}.fastq.gz`


## Usage on local system

    $ nextflow run main.nf -profile local

## Usage on an HPC with SLURM

To use the pipeline on a high performance cluster check if it supports the executor SLURM, the filesystem scratch2 and the modules Nextflow and Singularity (Ask your administrator for details). \
If requirements are met, create a batch file as given in the examples directory of this repository. Adjust parameters according to your cluster setup. Then run:

    $ sbatch example_batch_file.sh


## Pipeline parameters

### `-profile`

* This pipeline supports the two profiles `local` and `gwdg_cluster`
* `local` is using Docker as a container service.
* `slurm` is using Singularity as a container service and the executer SLURM (intended for HPC use).

Example:

    $ nextflow run main.nf -profile local

### `--reads`

* Specifies the location of the reads FASTQ file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* By default it is set to the endoMiner location: `$projectDir/examples/all_R{1,2}.fastq.gz`
* See above for naming convention of samples, replicates and pairs read files.

Example: 

    $ nextflow run main.nf --reads 'examples/all_R{1,2}.fastq.gz' -profile local

### `--endosymbiont_reference`

* The location of the endosymbiont reference genome fasta file.
* By default it is set to the endoMiner location: `$projectDir/examples/wNo_example_reference.fasta`.

Example:

    $ nextflow run main.nf --endosymbiont_reference 'examples/wNo_example_reference.fasta' -profile local

### `--output`

* Specifies desired output location for generated results.
* If the specified folder is not existant, it will be created.
* By default it is set to: `output`

Example:

    $ nextflow run main.nf --output 'output_folder/' -profile local

### `--max_memory`

* Specifies the amount of memory used by the pipeline in GB.
* By default it is set to: `128.GB`

Example:

    $ nextflow run main.nf --max_memory 128.GB -profile local

### `--max_cpus`

* Specifies the amount of CPU threads used by the pipeline.
* By default it is set to: `16`

Example:

    $ netflow run main.nf --max_cpus 16 -profile local

### `--max_time`

* Only relevant while executing on an HPC.
* Specifies the amount of time each process of the pipeline is allowed to run.
* Has to be in the format `1.h` for hours, `1.d` for days.
* By default it is set to: `1.d`

Example:

    $ nextflow run main.nf --max_time '1.d' -profile gwdg_cluster

### `--trim_*`

* This pipeline supports multiple TrimGalore! flags.
* `length`, `quality`, `adapter`, `phred64`, `clip_R1`, `three_prime_clip_R1`, `clip_R2` and `three_prime_clip_R2`
* To set them simply type `--trim_` followed by the desired flag.
* E.g. `--trim_length 55`
* By default `--trim_length` is set to `55` and `--trim_quality` to `20`. The other parameters are set to `None`.
* For details on TrimGalore! flags check the [TrimGalore! manual](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md).

Example:

    $ nextflow run main.nf --trim_length 55 -profile local

### Assembly flags `--kmers` and `--meta`

* This pipeline supports the SPAdes flags `--kmers` and `--meta`. 
* For details on `--meta` check the [SPAdes manual](https://github.com/ablab/spades#sec3.2). By default it is set to `true`.
* `--kmers` corresponds to the `--kmer` flag of SPAdes. For details check the [SPAdes manual](https://github.com/ablab/spades#k-mer-cardinality-estimating).
* To specify a single kmer length use for example `--kmers [111]`. For multiple kmer lenghts use for example `--kmers [71, 81, 91]`.

Example:

    $ nextflow run main.nf --meta true --kmers [71, 81, 91]


### `--min_blast_wordsize` and `--max_blast_wordsize`

* Sets the minimal and maximal word size for the mitogenome filtering process.
* HOSTMITOGENOMEFILTERING executes blastn with varying word sizes to find the mitogenome.
* If no single mitogenome is found after max_word_size is reached consider increasing it.
* By default `--min_blast_wordsize` is set to 11 and `--max_blast_wordsize` is set to 25.

Example:

    $ nextflow run main.nf --min_blast_wordsize 11 --max_blast_wordsize 25 -profile local

### `--help`

* Executes the help menu

Example:

    $ nextflow run main.nf --help

### `--version`

* Returns the version number of this pipeline

Example:

    $ nextflow run main.nf --version
