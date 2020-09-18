#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=regress_stdout.txt
#SBATCH --error=regress_stderr.txt
#SBATCH --time=1:12:00:00
#SBATCH --job-name=simulations
#SBATCH --mail-user=christopher.p.danko-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/danko/test
#
#################################################
module load R
Rscript tryboot-9.R > out.txt