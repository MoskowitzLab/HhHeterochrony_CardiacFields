
library(tidyverse)
library(slingshot)
library(tradeSeq)
library(BiocParallel)
library(stats)
library(parallel)
library(batchtools)

setwd("./2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_slingshot_condiments")
v6.1_sce <- readRDS("v6.1_sce_slingshot.rds")

#Number of knots: 9

# Set up parallel backend
param <- BiocParallel::bpparam()
param$workers <- 8

v6.1_sce_fitGAM <- fitGAM(v6.1_sce, conditions = factor(v6.1_sce$Genotype), nknots = 9, parallel = TRUE, BPPARAM = param)

saveRDS(v6.1_sce_fitGAM, "./2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq/v6.1_sce_fitGAM.rds")
