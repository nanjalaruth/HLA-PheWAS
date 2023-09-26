#!/bin/bash

#SBATCH --job-name='regression'
#SBATCH --output=reg-%j-stdout.log
#SBATCH --error=reg-%j-stderr.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=5

module load Nextflow/22.04.0

nextflow run main.nf -resume -profile slurm,singularity -c conf/test.config
