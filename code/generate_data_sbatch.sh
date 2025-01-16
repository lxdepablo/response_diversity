#!/bin/bash
#SBATCH -p atesting # Partition or queue
#SBATCH --job-name=generate_data # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=luis.depablo@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on 64 CPUs
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=0:10:00 # Time limit hrs:min:sec
#SBATCH --output=log_%j.out # Standard output and error log
#SBATCH --error=log_%j.err # %j inserts job number
pwd; hostname; date
echo "You've requested $SLURM_CPUS_ON_NODE core(s)."
date

rm *.err
rm *.out

module purge
module load R/4.4.0

Rscript generate_data.R
