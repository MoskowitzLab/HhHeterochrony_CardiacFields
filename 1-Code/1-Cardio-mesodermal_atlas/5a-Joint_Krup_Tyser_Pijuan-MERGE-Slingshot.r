#####################################################################################################################################################
################################                                                                                     ################################
################################ Krup-Pijuan-Tyser integrated Atlas - Merged dataset - Slingshot trajectory analysis ################################
################################                                                                                     ################################
#####################################################################################################################################################
# R/4.0.3
  library(dplyr)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(cowplot)
  library(sctransform)
  library(SeuratWrappers)
  library(viridis)
  library(slingshot)
  library(igraph)
  require(scales)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(grDevices)
  library(RColorBrewer)
  library(gam)
  library(Polychrome)
  library(circlize)
  library(ComplexHeatmap)
  library(DescTools)
  library(gridExtra)
  library(tidyHeatmap)
  library(msigdbr)
  library(fgsea)
  library(knitr)
  library(ggpubr)
  library(ggthemes)
  library(ggridges)
  library(ggrepel)
  library(forcats)
  library(factoextra)
  library(clusterProfiler)
  library(tidyverse)
  library(clusterExperiment)
  library(RSpectra)
  library(stringr)
  library(BBmisc)

#===============================================================================================================================================#
######################################################========      FHF trajectory analysis      ========########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E6.5-E7.5 & Hand1+ clusters

	# Time constraint
FHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"), invert=FALSE)
	# Marker gene expression constraint
FHF2_dataSub <- subset(x = FHF2_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,7,8,11,12), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/FHF2/UMAP_SubsetSelect_cells_FHF2.pdf", width=9, height=10)
print(DimPlot(FHF2_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/FHF2/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(FHF2_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c1
sce_FHF2_dataSub_start1_HARMONY30 <- as.SingleCellExperiment(FHF2_dataSub)
sce_FHF2_dataSub_start1_HARMONY30 <- slingshot(sce_FHF2_dataSub_start1_HARMONY30, clusterLabels = Idents(FHF2_dataSub), reducedDim = reducedDim(sce_FHF2_dataSub_start1_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="1", extend = 'n', stretch = 0)

saveRDS(sce_FHF2_dataSub_start1_HARMONY30, file = "./3-Merge3/Slingshot/FHF2/sce_FHF2_dataSub_SlingshotHARMONY30_STARTc1-extdN-Stretch0.RDS")
#sce_FHF2_dataSub_start1_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/FHF2/sce_FHF2_dataSub_SlingshotHARMONY30_STARTc1-extdN-Stretch0.RDS")

#Check the identified lineages
SlingshotDataSet(sce_FHF2_dataSub_start1_HARMONY30)
colnames(colData(sce_FHF2_dataSub_start1_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
pdf(file="./3-Merge3/Slingshot/FHF2/sce_FHF2_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc1-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(19,21)){
plotcol <- colors[cut(as.matrix(colData(sce_FHF2_dataSub_start1_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_FHF2_dataSub_start1_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_FHF2_dataSub_start1_HARMONY30)[u])))
}
dev.off()

#==# Selected cells only for pseudotime #==#
# Lineage 2 = FHF trajectory
GoodFHFCells <- rownames(na.omit(colData(sce_FHF2_dataSub_start1_HARMONY30)["slingPseudotime_2"]))

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0[["CellName"]] <- colnames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)
FHF_Slingshot_SelectedCells <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = CellName %in% GoodFHFCells, invert=FALSE)

meta <- na.omit(colData(sce_FHF2_dataSub_start1_HARMONY30)[,c("seurat_clusters","TimeStampMerged","slingPseudotime_2")])
colnames(meta) <- c("seurat_clusters", "TimeStampMerged", "dptval")
FHF_Slingshot_SelectedCells[["SlingshotPseudotime"]] <- meta$dptval[match(rownames(FHF_Slingshot_SelectedCells@meta.data), rownames(meta))]

### Save seurat object with only selected cells ###
saveRDS(FHF_Slingshot_SelectedCells, file = "./3-Merge3/Slingshot/FHF2/FHF_Slingshot_SelectedCells.RDS")
#FHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/FHF2/FHF_Slingshot_SelectedCells.RDS")
######

######################################
##  DEG along pseudotime - tradeSeq ##
######################################
meta <- na.omit(FHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Dataset")])
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Dataset")
meta_sort <- meta[order(meta$dptval),]

	# Choosing K
pseudotime <- as.matrix(meta_sort$dptval)
rownames(pseudotime) <- rownames(meta_sort)
colnames(pseudotime) <- "pseudotime"
cellWeights <- as.matrix(rep(1,length(pseudotime)))
rownames(cellWeights) <- rownames(meta_sort)
colnames(cellWeights) <- "cellWeights"
Dataset <- as.matrix(meta_sort$Dataset)
rownames(Dataset) <- rownames(meta_sort)
colnames(Dataset) <- "Dataset"

U <- model.matrix(~Dataset)
counts = as.matrix(FHF_Slingshot_SelectedCells@assays$RNA@counts)[FHF_Slingshot_SelectedCells@assays$SCT@var.features,rownames(meta_sort)] # only variable features
EvalK <- evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                   k=3:10,U = U, nGenes = 100, verbose = TRUE, plot = TRUE)

source("./0-R_codes/1-Functions/plot_evalutateK_results2.R")
pdf(file="./3-Merge3/Slingshot/FHF2/tradeSeq/EvalK.U-Dataset.FHF.pdf", width=10, height=7)
plot_evalutateK_results2(EvalK, k = NULL, aicDiff = 2)
dev.off()

GAM <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = TRUE)
				 
saveRDS(GAM, file = "./3-Merge3/Slingshot/FHF2/tradeSeq/GAM_FHF.RDS")
#GAM <- readRDS(file = "./3-Merge3/Slingshot/FHF2/tradeSeq/GAM_FHF.RDS")

assoRes <- associationTest(GAM) 
assoRes$padjusted <- p.adjust(assoRes$pvalue, method = "fdr")

FHF_DEGlist <- na.omit(assoRes[!(assoRes$padjusted > 0.05),])
FHF_DEGlist <- FHF_DEGlist[order(-FHF_DEGlist$meanLogFC), ]

write.csv(as.data.frame(FHF_DEGlist), file = "./3-Merge3/Slingshot/FHF2/tradeSeq/FHF_DEGlist.csv")

#===============================================================================================================================================#
######################################################========      JCF trajectory analysis      ========########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E7.25-E8.5 & Hand1+ clusters

	# Time constraint
JCF1_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
	# Marker gene expression constraint
JCF1_dataSub <- subset(x = JCF1_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,8,11), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/JCF1/UMAP_SubsetSelect_cells_JCF1.pdf", width=9, height=10)
print(DimPlot(JCF1_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/JCF1/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(JCF1_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c1
sce_JCF1_dataSub_start1_HARMONY30 <- as.SingleCellExperiment(JCF1_dataSub)
sce_JCF1_dataSub_start1_HARMONY30 <- slingshot(sce_JCF1_dataSub_start1_HARMONY30, clusterLabels = Idents(JCF1_dataSub), reducedDim = reducedDim(sce_JCF1_dataSub_start1_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="1", extend = 'n', stretch = 0)

saveRDS(sce_JCF1_dataSub_start1_HARMONY30, file = "./3-Merge3/Slingshot/JCF1/sce_JCF1_dataSub_SlingshotHARMONY30_STARTc1-extdN-Stretch0.RDS")
#sce_JCF1_dataSub_start1_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/JCF1/sce_JCF1_dataSub_SlingshotHARMONY30_STARTc1-extdN-Stretch0.RDS")

#Check the identified lineages
SlingshotDataSet(sce_JCF1_dataSub_start1_HARMONY30)
colnames(colData(sce_JCF1_dataSub_start1_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
pdf(file="./3-Merge3/Slingshot/JCF1/sce_JCF1_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc1-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(19,21)){
plotcol <- colors[cut(as.matrix(colData(sce_JCF1_dataSub_start1_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_JCF1_dataSub_start1_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_JCF1_dataSub_start1_HARMONY30)[u])))
}
dev.off()

#==# Selected cells only for pseudotime #==#
# Lineage 1 = JCF trajectory
GoodJCFCells <- rownames(na.omit(colData(sce_JCF1_dataSub_start1_HARMONY30)["slingPseudotime_1"]))

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0[["CellName"]] <- colnames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)
JCF_Slingshot_SelectedCells <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = CellName %in% GoodJCFCells, invert=FALSE)

meta <- na.omit(colData(sce_JCF1_dataSub_start1_HARMONY30)[,c("seurat_clusters","TimeStampMerged","slingPseudotime_1")])
colnames(meta) <- c("seurat_clusters", "TimeStampMerged", "dptval")
JCF_Slingshot_SelectedCells[["SlingshotPseudotime"]] <- meta$dptval[match(rownames(JCF_Slingshot_SelectedCells@meta.data), rownames(meta))]

### Save seurat object with only selected cells ###
saveRDS(JCF_Slingshot_SelectedCells, file = "./3-Merge3/Slingshot/JCF1/JCF_Slingshot_SelectedCells.RDS")
#JCF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/JCF1/JCF_Slingshot_SelectedCells.RDS")
######

######################################
##  DEG along pseudotime - tradeSeq ##
######################################
meta <- na.omit(JCF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Dataset")])
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Dataset")
meta_sort <- meta[order(meta$dptval),]

	# Choosing K
pseudotime <- as.matrix(meta_sort$dptval)
rownames(pseudotime) <- rownames(meta_sort)
colnames(pseudotime) <- "pseudotime"
cellWeights <- as.matrix(rep(1,length(pseudotime)))
rownames(cellWeights) <- rownames(meta_sort)
colnames(cellWeights) <- "cellWeights"
Dataset <- as.matrix(meta_sort$Dataset)
rownames(Dataset) <- rownames(meta_sort)
colnames(Dataset) <- "Dataset"

U <- model.matrix(~Dataset)
counts = as.matrix(JCF_Slingshot_SelectedCells@assays$RNA@counts)[JCF_Slingshot_SelectedCells@assays$SCT@var.features,rownames(meta_sort)] # only variable features
EvalK <- evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                   k=3:10,U = U, nGenes = 100, verbose = TRUE, plot = TRUE)

source("./0-R_codes/1-Functions/plot_evalutateK_results2.R")
pdf(file="./3-Merge3/Slingshot/JCF1/tradeSeq/EvalK.U-Dataset.JCF.pdf", width=10, height=7)
plot_evalutateK_results2(EvalK, k = NULL, aicDiff = 2)
dev.off()

GAM <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = TRUE)
				 
saveRDS(GAM, file = "./3-Merge3/Slingshot/JCF1/tradeSeq/GAM_JCF.RDS")
#GAM <- readRDS(file = "./3-Merge3/Slingshot/JCF1/tradeSeq/GAM_JCF.RDS")

assoRes <- associationTest(GAM) 
assoRes$padjusted <- p.adjust(assoRes$pvalue, method = "fdr")

JCF_DEGlist <- na.omit(assoRes[!(assoRes$padjusted > 0.05),])
JCF_DEGlist <- JCF_DEGlist[order(-JCF_DEGlist$meanLogFC), ]

write.csv(as.data.frame(JCF_DEGlist), file = "./3-Merge3/Slingshot/JCF1/tradeSeq/JCF_DEGlist.csv")

#===============================================================================================================================================#
######################################################========      SHF trajectory analysis      ========########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E7.25-E8.5 & Hand1+ clusters

	# Time constraint
SHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
	# Marker gene expression constraint
SHF2_dataSub <- subset(x = SHF2_dataPreSub, subset = seurat_clusters %in% c(1,6,9), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/SHF2/UMAP_SubsetSelect_cells_SHF2.pdf", width=9, height=10)
print(DimPlot(SHF2_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/SHF2/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(SHF2_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c1
sce_SHF2_dataSub_start1_HARMONY30 <- as.SingleCellExperiment(SHF2_dataSub)
sce_SHF2_dataSub_start1_HARMONY30 <- slingshot(sce_SHF2_dataSub_start1_HARMONY30, clusterLabels = Idents(SHF2_dataSub), reducedDim = reducedDim(sce_SHF2_dataSub_start1_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="1", end.clus = c("6","9"), extend = 'n', stretch = 0)

saveRDS(sce_SHF2_dataSub_start1_HARMONY30, file = "./3-Merge3/Slingshot/SHF2/sce_SHF2_dataSub_SlingshotHARMONY30_STARTc1-ENDc6-9-extdN-Stretch0.RDS")
#sce_SHF2_dataSub_start1_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/SHF2/sce_SHF2_dataSub_SlingshotHARMONY30_STARTc1-ENDc6-9-extdN-Stretch0.RDS")

#Check the identified lineages and Curves
SlingshotDataSet(sce_SHF2_dataSub_start1_HARMONY30)
colnames(colData(sce_SHF2_dataSub_start1_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))
pdf(file="./3-Merge3/Slingshot/SHF2/sce_SHF2_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc1-ENDc6-9-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(19,20)){
plotcol <- colors[cut(as.matrix(colData(sce_SHF2_dataSub_start1_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_SHF2_dataSub_start1_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_SHF2_dataSub_start1_HARMONY30)[u])))
#lines(SlingshotDataSet(sce_SHF2_dataSub_start1_HARMONY30), lwd=2, col='black')
}
dev.off()

#==# Selected cells only for pseudotime #==#
# Lineage 1 = aSHF trajectory
GoodaSHFCells <- rownames(na.omit(colData(sce_SHF2_dataSub_start1_HARMONY30)["slingPseudotime_1"]))

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0[["CellName"]] <- colnames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)
aSHF_Slingshot_SelectedCells <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = CellName %in% GoodaSHFCells, invert=FALSE)

meta <- na.omit(colData(sce_SHF2_dataSub_start1_HARMONY30)[,c("seurat_clusters","TimeStampMerged","slingPseudotime_1")])
colnames(meta) <- c("seurat_clusters", "TimeStampMerged", "dptval")
aSHF_Slingshot_SelectedCells[["SlingshotPseudotime"]] <- meta$dptval[match(rownames(aSHF_Slingshot_SelectedCells@meta.data), rownames(meta))]

### Save seurat object with only selected cells ###
saveRDS(aSHF_Slingshot_SelectedCells, file = "./3-Merge3/Slingshot/SHF2/aSHF_Slingshot_SelectedCells.RDS")
#aSHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/SHF2/aSHF_Slingshot_SelectedCells.RDS")
######

#==# Selected cells only for pseudotime #==#
# Lineage 2 = pSHF trajectory
GoodpSHFCells <- rownames(na.omit(colData(sce_SHF2_dataSub_start1_HARMONY30)["slingPseudotime_2"]))

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0[["CellName"]] <- colnames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)
pSHF_Slingshot_SelectedCells <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = CellName %in% GoodpSHFCells, invert=FALSE)

meta <- na.omit(colData(sce_SHF2_dataSub_start1_HARMONY30)[,c("seurat_clusters","TimeStampMerged","slingPseudotime_2")])
colnames(meta) <- c("seurat_clusters", "TimeStampMerged", "dptval")
pSHF_Slingshot_SelectedCells[["SlingshotPseudotime"]] <- meta$dptval[match(rownames(pSHF_Slingshot_SelectedCells@meta.data), rownames(meta))]

### Save seurat object with only selected cells ###
saveRDS(pSHF_Slingshot_SelectedCells, file = "./3-Merge3/Slingshot/SHF2/pSHF_Slingshot_SelectedCells.RDS")
#pSHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/SHF2/pSHF_Slingshot_SelectedCells.RDS")
######

######################################
##  DEG along pseudotime - tradeSeq ##
######################################
################ aSHF ################
meta <- na.omit(aSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Dataset")])
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Dataset")
meta_sort <- meta[order(meta$dptval),]

	# Choosing K
pseudotime <- as.matrix(meta_sort$dptval)
rownames(pseudotime) <- rownames(meta_sort)
colnames(pseudotime) <- "pseudotime"
cellWeights <- as.matrix(rep(1,length(pseudotime)))
rownames(cellWeights) <- rownames(meta_sort)
colnames(cellWeights) <- "cellWeights"
Dataset <- as.matrix(meta_sort$Dataset)
rownames(Dataset) <- rownames(meta_sort)
colnames(Dataset) <- "Dataset"

U <- model.matrix(~Dataset)
counts = as.matrix(aSHF_Slingshot_SelectedCells@assays$RNA@counts)[aSHF_Slingshot_SelectedCells@assays$SCT@var.features,rownames(meta_sort)] # only variable features
EvalK <- evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                   k=3:10,U = U, nGenes = 100, verbose = TRUE, plot = TRUE)

source("./0-R_codes/1-Functions/plot_evalutateK_results2.R")
pdf(file="./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/EvalK.U-Dataset.aSHF.pdf", width=10, height=7)
plot_evalutateK_results2(EvalK, k = NULL, aicDiff = 2)
dev.off()

GAM <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = TRUE)
				 
saveRDS(GAM, file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/GAM_aSHF.RDS")
#GAM <- readRDS(file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/GAM_aSHF.RDS")

assoRes <- associationTest(GAM) 
assoRes$padjusted <- p.adjust(assoRes$pvalue, method = "fdr")

aSHF_DEGlist <- na.omit(assoRes[!(assoRes$padjusted > 0.05),])
aSHF_DEGlist <- aSHF_DEGlist[order(-aSHF_DEGlist$meanLogFC), ]

write.csv(as.data.frame(aSHF_DEGlist), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/aSHF_DEGlist.csv")

################ pSHF ################
meta <- na.omit(pSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Dataset")])
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Dataset")
meta_sort <- meta[order(meta$dptval),]

	# Choosing K
pseudotime <- as.matrix(meta_sort$dptval)
rownames(pseudotime) <- rownames(meta_sort)
colnames(pseudotime) <- "pseudotime"
cellWeights <- as.matrix(rep(1,length(pseudotime)))
rownames(cellWeights) <- rownames(meta_sort)
colnames(cellWeights) <- "cellWeights"
Dataset <- as.matrix(meta_sort$Dataset)
rownames(Dataset) <- rownames(meta_sort)
colnames(Dataset) <- "Dataset"

U <- model.matrix(~Dataset)
counts = as.matrix(pSHF_Slingshot_SelectedCells@assays$RNA@counts)[pSHF_Slingshot_SelectedCells@assays$SCT@var.features,rownames(meta_sort)] # only variable features
EvalK <- evaluateK(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                   k=3:10,U = U, nGenes = 100, verbose = TRUE, plot = TRUE)

source("./0-R_codes/1-Functions/plot_evalutateK_results2.R")
pdf(file="./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/EvalK.U-Dataset.pSHF.pdf", width=10, height=7)
plot_evalutateK_results2(EvalK, k = NULL, aicDiff = 2)
dev.off()

GAM <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = TRUE)
				 
saveRDS(GAM, file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/GAM_pSHF.RDS")
#GAM <- readRDS(file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/GAM_pSHF.RDS")

assoRes <- associationTest(GAM) 
assoRes$padjusted <- p.adjust(assoRes$pvalue, method = "fdr")

pSHF_DEGlist <- na.omit(assoRes[!(assoRes$padjusted > 0.05),])
pSHF_DEGlist <- pSHF_DEGlist[order(-pSHF_DEGlist$meanLogFC), ]

write.csv(as.data.frame(pSHF_DEGlist), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/pSHF_DEGlist.csv")

#===============================================================================================================================================#
###########################################========      Extended FHF trajectory analysis      ========##########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E6.5-E7.5 & Hand1+ clusters
# Add 'early' cell types/clusters

ExtendedFHF_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"), invert=FALSE)
ExtendedFHF_dataSub <- subset(x = ExtendedFHF_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,7,8,11,12,10,3,16,17,15,0), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/ExtendedFHF/UMAP_SubsetSelect_cells_ExtendedFHF.pdf", width=9, height=10)
print(DimPlot(ExtendedFHF_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/ExtendedFHF/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(ExtendedFHF_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c0
sce_ExtendedFHF_dataSub_start0_HARMONY30 <- as.SingleCellExperiment(ExtendedFHF_dataSub)
sce_ExtendedFHF_dataSub_start0_HARMONY30 <- slingshot(sce_ExtendedFHF_dataSub_start0_HARMONY30, clusterLabels = Idents(ExtendedFHF_dataSub), reducedDim = reducedDim(sce_ExtendedFHF_dataSub_start0_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="0", extend = 'n', stretch = 0)

saveRDS(sce_ExtendedFHF_dataSub_start0_HARMONY30, file = "./3-Merge3/Slingshot/ExtendedFHF/sce_ExtendedFHF_dataSub_SlingshotHARMONY30_STARTc0-extdN-Stretch0.RDS")
#sce_ExtendedFHF_dataSub_start0_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/ExtendedFHF/sce_ExtendedFHF_dataSub_SlingshotHARMONY30_STARTc0-extdN-Stretch0.RDS")

#Check the identified lineages and Curves
SlingshotDataSet(sce_ExtendedFHF_dataSub_start0_HARMONY30)
colnames(colData(sce_ExtendedFHF_dataSub_start0_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))

pdf(file="./3-Merge3/Slingshot/ExtendedFHF/sce_ExtendedFHF_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc0-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(20,26)){
plotcol <- colors[cut(as.matrix(colData(sce_ExtendedFHF_dataSub_start0_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_ExtendedFHF_dataSub_start0_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_ExtendedFHF_dataSub_start0_HARMONY30)[u])))
#lines(SlingshotDataSet(sce_ExtendedFHF_dataSub_start0_HARMONY30), lwd=2, col='black')
}
dev.off()

#===============================================================================================================================================#
###########################################========      Extended JCF trajectory analysis      ========##########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E7.25-E8.5 & Hand1+ clusters
# Add 'early' cell types/clusters

ExtendedJCF_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
ExtendedJCF_dataSub <- subset(x = ExtendedJCF_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,8,11,10,3,16,17,15,0), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/ExtendedJCF/UMAP_SubsetSelect_cells_ExtendedJCF.pdf", width=9, height=10)
print(DimPlot(ExtendedJCF_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/ExtendedJCF/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(ExtendedJCF_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c0
sce_ExtendedJCF_dataSub_start0_HARMONY30 <- as.SingleCellExperiment(ExtendedJCF_dataSub)
sce_ExtendedJCF_dataSub_start0_HARMONY30 <- slingshot(sce_ExtendedJCF_dataSub_start0_HARMONY30, clusterLabels = Idents(ExtendedJCF_dataSub), reducedDim = reducedDim(sce_ExtendedJCF_dataSub_start0_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="0", extend = 'n', stretch = 0)

saveRDS(sce_ExtendedJCF_dataSub_start0_HARMONY30, file = "./3-Merge3/Slingshot/ExtendedJCF/sce_ExtendedJCF_dataSub_SlingshotHARMONY30_STARTc0-extdN-Stretch0.RDS")
#sce_ExtendedJCF_dataSub_start0_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/ExtendedJCF/sce_ExtendedJCF_dataSub_SlingshotHARMONY30_STARTc0-extdN-Stretch0.RDS")

#Check the identified lineages and Curves
SlingshotDataSet(sce_ExtendedJCF_dataSub_start0_HARMONY30)
colnames(colData(sce_ExtendedJCF_dataSub_start0_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))

pdf(file="./3-Merge3/Slingshot/ExtendedJCF/sce_ExtendedJCF_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc0-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(20,24)){
plotcol <- colors[cut(as.matrix(colData(sce_ExtendedJCF_dataSub_start0_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_ExtendedJCF_dataSub_start0_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_ExtendedJCF_dataSub_start0_HARMONY30)[u])))
#lines(SlingshotDataSet(sce_ExtendedJCF_dataSub_start0_HARMONY30), lwd=2, col='black')
}
dev.off()

#===============================================================================================================================================#
###########################################========      Extended SHF trajectory analysis      ========##########################################
#===============================================================================================================================================#
######################################
## Slingshot trajectory calculation ##
######################################
# Constraints: E7.25-E8.5 & Isl1+ clusters
# Add 'early' cell types/clusters

ExtendedSHF_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
ExtendedSHF_dataSub <- subset(x = ExtendedSHF_dataPreSub, subset = seurat_clusters %in% c(1,6,9,10,3,16,17,15,0), invert=FALSE)

pdf(file="./3-Merge3/Slingshot/ExtendedSHF/UMAP_SubsetSelect_cells_ExtendedSHF.pdf", width=9, height=10)
print(DimPlot(ExtendedSHF_dataSub, reduction = "umap", label = T, raster = FALSE))
dev.off()

pdf(file="./3-Merge3/Slingshot/ExtendedSHF/ElbowPlot_Harmony.pdf",width=5,height=5)
ElbowPlot(ExtendedSHF_dataSub, ndims = 50, reduction = "harmony")
dev.off()
# 30 components

##START c0
sce_ExtendedSHF_dataSub_start1_HARMONY30 <- as.SingleCellExperiment(ExtendedSHF_dataSub)
sce_ExtendedSHF_dataSub_start1_HARMONY30 <- slingshot(sce_ExtendedSHF_dataSub_start1_HARMONY30, clusterLabels = Idents(ExtendedSHF_dataSub), reducedDim = reducedDim(sce_ExtendedSHF_dataSub_start1_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="0", end.clus = c("6","9"), extend = 'n', stretch = 0)

saveRDS(sce_ExtendedSHF_dataSub_start1_HARMONY30, file = "./3-Merge3/Slingshot/ExtendedSHF/sce_ExtendedSHF_dataSub_SlingshotHARMONY30_STARTc0-ENDc6-9-extdN-Stretch0.RDS")
#sce_ExtendedSHF_dataSub_start1_HARMONY30 <- readRDS(file = "./3-Merge3/Slingshot/ExtendedSHF/sce_ExtendedSHF_dataSub_SlingshotHARMONY30_STARTc0-ENDc6-9-extdN-Stretch0.RDS")

#Check the identified lineages and Curves
SlingshotDataSet(sce_ExtendedSHF_dataSub_start1_HARMONY30)

colnames(colData(sce_ExtendedSHF_dataSub_start1_HARMONY30))

	# Plotting pseudotimes
colors <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100))

pdf(file="./3-Merge3/Slingshot/ExtendedSHF/sce_ExtendedSHF_dataSub_SlingshotHARMONY30_One-functionMethod_Pseudotimes_HARMONY30_STARTc0-ENDc6-9-extdN-Stretch0nolines.pdf", width=12, height=12)
for(u in seq(20,24)){
plotcol <- colors[cut(as.matrix(colData(sce_ExtendedSHF_dataSub_start1_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_ExtendedSHF_dataSub_start1_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_ExtendedSHF_dataSub_start1_HARMONY30)[u])))
#lines(SlingshotDataSet(sce_ExtendedSHF_dataSub_start1_HARMONY30), lwd=2, col='black')
}
dev.off()
