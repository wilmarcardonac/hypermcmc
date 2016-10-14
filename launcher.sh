#!/bin/sh
#SBATCH --cpus-per-task=4
#SBATCH --job-name=hyper-MCMC
#SBATCH --ntasks=1
#SBATCH --time=0-12:00:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
##SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out

srun ./mcmc