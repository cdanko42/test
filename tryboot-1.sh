#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=1024
#SBATCH --output=regress_stdout.txt
#SBATCH --error=regress_stderr.txt
<<<<<<< HEAD
#SBATCH --time=1:12:00:00
=======
#SBATCH --time=15:00:00
>>>>>>> 011ee7261d33b9a1d73edbb31f8c50ea8dbb80a1
#SBATCH --job-name=simulations
#SBATCH --mail-user=christopher.p.danko-1@ou.edu
#SBATCH --mail-type=ALL
#SBATCH --chdir=/home/danko/test
#
#################################################
module load R
<<<<<<< HEAD
Rscript tryboot-1.R > out.txt
=======
Rscript tryboot.R > out.txt
>>>>>>> 011ee7261d33b9a1d73edbb31f8c50ea8dbb80a1
