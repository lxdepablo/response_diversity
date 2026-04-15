#!/bin/bash
#SBATCH -p amilan # Partition or queue
#SBATCH --qos=normal
#SBATCH --job-name=generate_data # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
<<<<<<< HEAD
#SBATCH --ntasks=16 # Run on 64 CPUs
#SBATCH --mem=12gb # Memory limit
#SBATCH --time=2:00:00 # Time limit hrs:min:sec
=======
#SBATCH --ntasks=64 # Run on 64 CPUs
#SBATCH --mem=64gb # Memory limit
#SBATCH --time=1:00:00 # Time limit hrs:min:sec
>>>>>>> 668e8cb5b0b315e584f9206ee57f97c3f8472842
#SBATCH --output=log.out # Standard output and error log
#SBATCH --error=log.err # %j inserts job number

pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

module purge
module load miniforge
mamba activate tidyverse_env

Rscript generate_data.R
