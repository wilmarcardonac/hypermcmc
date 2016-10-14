#!/bin/sh
#SBATCH --cpus-per-task=4
#SBATCH --job-name=R11-LMC
#SBATCH --ntasks=1
#SBATCH --time=0-8:00:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out

srun ./mcmc