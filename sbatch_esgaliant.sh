#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=esgaliant_env
#SBATCH --mem-per-cpu=1GB
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>
#SBATCH -p defq
#SBATCH --nodelist=esplhpc-cp047

module load nextflow
module load apptainer
nextflow run main.nf -c configs/run.config