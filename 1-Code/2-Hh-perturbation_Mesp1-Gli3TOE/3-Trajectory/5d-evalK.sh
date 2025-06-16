#!/bin/bash

#SBATCH --job-name=evalk
#SBATCH --account=pi-imoskowitz
#SBATCH --output=evalk.out
#SBATCH --error=evalk.err
#SBATCH --time=20:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wmikolajczak@rcc.uchicago.edu
cd /project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq

module load 
module load gcc 
module load ImageMagick 
module load R/4.3.1
cd /project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq/sbatch
Rscript evalK.R
