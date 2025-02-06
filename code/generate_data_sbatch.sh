#!/bin/bash
#SBATCH -p amilan # Partition or queue
#SBATCH --job-name=generate_data # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=64 # Run on 64 CPUs
#SBATCH --mem=64gb # Memory limit
#SBATCH --time=2:00:00 # Time limit hrs:min:sec
#SBATCH --output=log.out # Standard output and error log
#SBATCH --error=log.err # %j inserts job number

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

module purge
module load mambaforge
mamba activate tidyverse_env

Rscript generate_data.R
