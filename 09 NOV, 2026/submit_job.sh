#!/bin/bash
#SBATCH -J large_R_program
#SBATCH -A genacc_q
#SBATCH -p genacc_q
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 01-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=st24g@fsu.edu

module purge
module load R

export R_LIBS_USER=/gpfs/home/st24g/R/x86_64-redhat-linux-gnu-library/4.5

Rscript "simulation.R"
