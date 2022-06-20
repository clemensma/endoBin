#!/bin/bash
#SBATCH --job-name=your_job_name
#SBATCH --partition=medium
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --constraint=scratch2
#SBATCH --time=5-00:00:00
#SBATCH --mem=32G
#SBATCH --output=/path/to/log_file.txt
#SBATCH --error=/path/to/error_file.txt

module purge
module load nextflow
module load singularity

nextflow run '/path/to/endoMiner/main.nf'\
  --reads '/path/to/endoMiner/examples/all_R{1,2}\.fastq\.gz'\
  --endosymbiont_reference '/path/to/endoMiner/examples/wNo_example_reference.fasta'\
  -profile slurm\
  --memory 1024\
  --threads 32\
  --time '47.h'\
  --output '/path/to/output/folder'
