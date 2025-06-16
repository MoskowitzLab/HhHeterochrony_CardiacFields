#!/bin/bash

#SBATCH --job-name=fitgam
#SBATCH --account=pi-imoskowitz
#SBATCH --output=fitgam.out
#SBATCH --error=fitgam.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=700G
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wmikolajczak@rcc.uchicago.edu
cd /project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq

module load 
module load gcc 
module load ImageMagick 
module load R/4.3.1
cd sbatch

echo "Initial Memory Usage:"
free -h

Rscript fitgam.R

echo "Final Memory Usage:"
free -h

