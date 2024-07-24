#!/bin/bash
#SBATCH -p short # Partition or queue
#SBATCH --job-name=generate_data # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=64 # Run on 16 CPUs
#SBATCH --mem=50gb # Memory limit
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/lude8513/log_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/lude8513/log_%j.err # %j inserts job number
pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

module purge
module load R/4.3.1

Rscript generate_data.R