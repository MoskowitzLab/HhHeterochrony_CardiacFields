library(slingshot)
library(tradeSeq)
library(tidyverse)

setwd("/project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_slingshot_condiments")

v6.1_sce <- readRDS("v6.1_sce_slingshot.rds")
pdf("/project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq/EvalKoutput.pdf")
icMat <- tradeSeq::evaluateK(v6.1_sce, conditions= as.factor(v6.1_sce$Genotype), k=3:15)
dev.off()

saveRDS(icMat, "/project/imoskowitz/wiktor/2-projects/2-Mesp1creGli3T_scRNAseq_Guzetta/V6/v6.1_tradeseq/icMat.rds")

