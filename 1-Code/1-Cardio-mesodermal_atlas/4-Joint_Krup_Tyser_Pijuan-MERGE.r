#####################################################################################################################################################
################################################                                                     ################################################
################################################ Krup-Pijuan-Tyser integrated Atlas - Merged dataset ################################################
################################################                                                     ################################################
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
######################################################========      Merge 1      ========#######################################################
#===============================================================================================================================================#
## Data loading ##
SeuratTyser <- readRDS(file = "./3-PublicDatasets/8-TyserJCF_2021/Integration/SeuratTyser_SCT.RDS")
v1.2.6.5_SUBSET1_harmony_res1.0 <- readRDS(file = "./Krup_GSE208153_scRNAseq/0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")
SeuratGastrulationMeso <- readRDS(file = "./MouseGastrulationData/SeuratGastrulationMeso_SCT.RDS")

######################################
## Data cleaning and homogenization ##
######################################
	# Tyser
SeuratTyser@meta.data$Batch <- SeuratTyser@meta.data$batch
SeuratTyser@meta.data$Dataset <- "Tyser.2021"
SeuratTyser@meta.data$Timepoint <- SeuratTyser@meta.data$stage
SeuratTyser@meta.data$Old.groups <- SeuratTyser@meta.data$clusterAnn

SeuratTyser@meta.data$batch <- NULL
SeuratTyser@meta.data$sampleID <- NULL
SeuratTyser@meta.data$date <- NULL
SeuratTyser@meta.data$plate <- NULL
SeuratTyser@meta.data$well <- NULL
SeuratTyser@meta.data$stage <- NULL
SeuratTyser@meta.data$side <- NULL
SeuratTyser@meta.data$ENA <- NULL
SeuratTyser@meta.data$cluster <- NULL
SeuratTyser@meta.data$clusterAnn <- NULL
SeuratTyser@meta.data$clusterCol <- NULL
SeuratTyser@meta.data$old.ident <- NULL

head(SeuratTyser@meta.data)

	#Krup
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Batch <- v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Library.ident
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Dataset <- "Krup.2023"
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Timepoint <- v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$gem.group
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Old.groups <- v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$seurat_clusters

v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$gem.group <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$percent.mt <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$old.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$Library.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$SCT_snn_res.1 <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$seurat_clusters <- NULL

head(v1.2.6.5_SUBSET1_harmony_res1.0@meta.data)

	# Pijuan
SeuratGastrulationMeso@meta.data$Batch <- SeuratGastrulationMeso@meta.data$sequencing.batch
SeuratGastrulationMeso@meta.data$Dataset <- "Pijuan.2019"
SeuratGastrulationMeso@meta.data$Timepoint <- SeuratGastrulationMeso@meta.data$stage
SeuratGastrulationMeso@meta.data$Old.groups <- SeuratGastrulationMeso@meta.data$celltype

SeuratGastrulationMeso@meta.data$cell <- NULL
SeuratGastrulationMeso@meta.data$barcode <- NULL
SeuratGastrulationMeso@meta.data$sample <- NULL
SeuratGastrulationMeso@meta.data$pool <- NULL
SeuratGastrulationMeso@meta.data$stage <- NULL
SeuratGastrulationMeso@meta.data$sequencing.batch <- NULL
SeuratGastrulationMeso@meta.data$theiler <- NULL
SeuratGastrulationMeso@meta.data$doub.density <- NULL
SeuratGastrulationMeso@meta.data$doublet <- NULL
SeuratGastrulationMeso@meta.data$cluster <- NULL
SeuratGastrulationMeso@meta.data$cluster.sub <- NULL
SeuratGastrulationMeso@meta.data$cluster.stage <- NULL
SeuratGastrulationMeso@meta.data$cluster.theiler <- NULL
SeuratGastrulationMeso@meta.data$stripped <- NULL
SeuratGastrulationMeso@meta.data$celltype <- NULL
SeuratGastrulationMeso@meta.data$colour <- NULL
SeuratGastrulationMeso@meta.data$sizeFactor <- NULL
SeuratGastrulationMeso@meta.data$Limited_Hh_score_1 <- NULL
SeuratGastrulationMeso@meta.data$Me4_1 <- NULL
SeuratGastrulationMeso@meta.data$Me7_1 <- NULL
SeuratGastrulationMeso@meta.data$Me5_1 <- NULL
SeuratGastrulationMeso@meta.data$Me6_1 <- NULL
SeuratGastrulationMeso@meta.data$old.ident <- NULL

head(SeuratGastrulationMeso@meta.data)

############################
## First quality controls ##
############################
Vln_Krup <- VlnPlot(v1.2.6.5_SUBSET1_harmony_res1.0, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Tyser <- VlnPlot(SeuratTyser, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Pijuan <- VlnPlot(SeuratGastrulationMeso, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")

plot_Krup1 <- FeatureScatter(v1.2.6.5_SUBSET1_harmony_res1.0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Tyser1 <- FeatureScatter(SeuratTyser, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Pijuan1 <- FeatureScatter(SeuratGastrulationMeso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(file="./1-FirstSeparate-QC_Plots_raw.pdf", width=9, height=10)
print(Vln_Krup)
print(Vln_Tyser)
print(Vln_Pijuan)
print(plot_Krup1)
print(plot_Tyser1)
print(plot_Pijuan1)
dev.off()

####################
## Merge datasets ##
####################

Merge_Krup_Tyser_Pijuan <- merge(v1.2.6.5_SUBSET1_harmony_res1.0,c(SeuratTyser,SeuratGastrulationMeso))

VarFeatures <- rbind(data.frame(x=VariableFeatures(v1.2.6.5_SUBSET1_harmony_res1.0)),data.frame(x=VariableFeatures(SeuratTyser)),data.frame(x=VariableFeatures(SeuratGastrulationMeso)))
VarFeatures_List <- unique(VarFeatures$x)

VariableFeatures(Merge_Krup_Tyser_Pijuan) <- VarFeatures_List

# 49,668 cells - 26,845 SCT features / 6,104 variable features

saveRDS(Merge_Krup_Tyser_Pijuan, file = "./Merge_Krup_Tyser_Pijuan.RDS")
#Merge_Krup_Tyser_Pijuan <- readRDS(file = "./Merge_Krup_Tyser_Pijuan.RDS")

############################
## Merge quality controls ##
############################
Vln_Merge <- VlnPlot(Merge_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, group.by = "Dataset")
Vln_MergeSCT <- VlnPlot(Merge_Krup_Tyser_Pijuan, features = c("nFeature_SCT", "nCount_SCT"), ncol = 3, group.by = "Dataset")
plot_Merge <- FeatureScatter(Merge_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_MergeSCT <- FeatureScatter(Merge_Krup_Tyser_Pijuan, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")

pdf(file="./1-Merge1/1-QC_Plots_Merge.pdf", width=9, height=10)
print(Vln_Merge)
print(Vln_MergeSCT)
print(plot_Merge)
print(plot_MergeSCT)
dev.off()

###################################################
## PCA, Harmony Integration, UMAP and clustering ##
###################################################
Merge_Krup_Tyser_Pijuan <- RunPCA(Merge_Krup_Tyser_Pijuan, verbose = TRUE)

pdf(file="./1-Merge1/Merge_Krup_Tyser_Pijuan_ElbowPlot-PCA.pdf", width=9, height=10)
ElbowPlot(Merge_Krup_Tyser_Pijuan, ndims = 50, reduction = "pca")
dev.off()
# 30 components

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
Merge_Krup_Tyser_Pijuan_harmony2 <- RunHarmony(Merge_Krup_Tyser_Pijuan, 
				group.by.vars = c("Batch","Dataset"), 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

print("Integration / Batch effect correction done!")

# UMAP calculation and clustering
Merge_Krup_Tyser_Pijuan_harmony2 <- RunUMAP(Merge_Krup_Tyser_Pijuan_harmony2, dims = 1:30, reduction = "harmony", verbose = TRUE)
Merge_Krup_Tyser_Pijuan_harmony2 <- FindNeighbors(Merge_Krup_Tyser_Pijuan_harmony2,reduction = "harmony", dims = 1:30, verbose = TRUE)
Merge_Krup_Tyser_Pijuan_harmony2_res1.0 <- FindClusters(Merge_Krup_Tyser_Pijuan_harmony2, verbose = TRUE, resolution = 1.0)

pdf(file="./1-Merge1/UMAP_Merge_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T, raster = FALSE))
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Dataset", raster = FALSE))
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Batch", raster = FALSE))
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "orig.ident", raster = FALSE)+ NoLegend())
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Timepoint", raster = FALSE))
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE))
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE)+ NoLegend())
dev.off()

# Add Log count per million normalization (RNA slot) for gene expression visualization and analysis
Merge_Krup_Tyser_Pijuan_harmony2_res1.0 <- NormalizeData(
  Merge_Krup_Tyser_Pijuan_harmony2_res1.0,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1000000,
  verbose = TRUE
)

# Homogeneate timestamps
Merge_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(TimeStamp = case_when(
    Timepoint %in% "E60" ~ "E6.0",
	Timepoint %in% "lateE60" ~ "lateE6.0",
	Timepoint %in% "E65" ~ "E6.5",
	Timepoint %in% "lateE65" ~ "lateE6.5",
	Timepoint %in% "E70" ~ "E7.0",
	Timepoint %in% "earlyE75" ~ "earlyE7.5",
	Timepoint %in% "lateE75" ~ "lateE7.5",
	Timepoint %in% "E775" ~ "E7.75",
	TRUE ~ Timepoint
  ))

# Time Stamp UMAP
pdf(file="./1-Merge1/TimeStamp-UMAP_Merge_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
for(Day in unique(Merge_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimeStamp)){
cells <- WhichCells(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimeStamp ==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", group.by = "TimeStamp", cells.highlight=Cells, cols.highlight=c("#4472C4"),sizes.highlight = 0.6, pt.size = 0.4, cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

print("Time Stamp UMAP done!")

#######################################################################
## Gene expression analysis - need to switch from SCT to RNA for CPM ##
#######################################################################
TmpObj_GeneExpr <- Merge_Krup_Tyser_Pijuan_harmony2_res1.0
DefaultAssay(TmpObj_GeneExpr) <- "RNA"

# Dot plot and Feature plots QCs
gene_list <- c("Fgf5","Isl1","Osr1","Tbx1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Asb2","Nkx2-5","Tbx5","Gata4","Gata5","Gata6","Gli1","Foxf1","Hhip","Ptch1","Fgf8","Fgf10","Mab21l2","Hand1","Smarcd3","Ai14-transgene","eGFP-transgene","Myl2","Hcn4","Hopx")
pdf(file = "./1-Merge1/MarkerGenes.v1_Merge_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(TmpObj_GeneExpr, features = i) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

# Version data save
saveRDS(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, file = "./1-Merge1/Merge_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
#Merge_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./1-Merge1/Merge_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

print("Version saved!")

################################
## Top markers identification ##
################################
Merge_Krup_Tyser_Pijuan_harmony2_res1.0 <- PrepSCTFindMarkers(Merge_Krup_Tyser_Pijuan_harmony2_res1.0)
Merge_Krup_Tyser_Pijuan_harmony2_res1.0.markers <- FindAllMarkers(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6,recorrect_umi = FALSE)
write.csv(Merge_Krup_Tyser_Pijuan_harmony2_res1.0.markers, file = "./1-Merge1/Merge_Krup_Tyser_Pijuan_harmony2_res1.0_all-markers.csv")

# Marker heatmap
Merge_Krup_Tyser_Pijuan_harmony2_res1.0.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = "./1-Merge1/TopMarkers_Merge_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 11, height = 11)
print(DoHeatmap(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, features = top15$gene) + NoLegend())
dev.off()

pdf(file = "./1-Merge1/DotPlot_Merge_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 15, height = 40)
print(DotPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

#===============================================================================================================================================#
######################################################========      Subset 1      ========#######################################################
#===============================================================================================================================================#
# SUBSET preparation - need to remove clusters: 9,10,11,17,19,20,21,22,24,25,26 (from Merge_Krup_Tyser_Pijuan_harmony2_res1.0)
# Need to obtain cell names to remove them from individual dataset followed by separate SCT before merging/Harmony integration

CellsToRemove <- as.data.frame(Merge_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data[Merge_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$seurat_clusters %in% c(9,10,11,17,19,20,21,22,24,25,26),c("seurat_clusters","Dataset")])

CellsToRemove_Krup <- rownames(CellsToRemove[CellsToRemove$Dataset == "Krup.2023",])
write.csv(CellsToRemove_Krup, file = "./1-Merge1/PrepSub1_CellsToRemove_Krup.csv")

CellsToRemove_Tyser <- rownames(CellsToRemove[CellsToRemove$Dataset == "Tyser.2021",])
write.csv(CellsToRemove_Tyser, file = "./1-Merge1/PrepSub1_CellsToRemove_Tyser.csv")

CellsToRemove_Pijuan <- rownames(CellsToRemove[CellsToRemove$Dataset == "Pijuan.2019",])
write.csv(CellsToRemove_Pijuan, file = "./1-Merge1/PrepSub1_CellsToRemove_Pijuan.csv")

#===============================================================================================================================================#
######################################################========      Merge 2      ========#######################################################
#===============================================================================================================================================#
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./Krup_GSE208153_scRNAseq/0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan.RDS")
SeuratTyser_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./3-PublicDatasets/8-TyserJCF_2021/Integration/SeuratTyser_subForMerge_Krup_Tyser_Pijuan.RDS")
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./MouseGastrulationData/SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan.RDS")

######################################
## Data cleaning and homogenization ##
######################################
	# Tyser
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$Batch <- SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$batch
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$Dataset <- "Tyser.2021"
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$Timepoint <- SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$stage
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$Old.groups <- SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$clusterAnn

SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$batch <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$sampleID <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$date <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$plate <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$well <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$stage <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$side <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$ENA <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$cluster <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$clusterAnn <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$clusterCol <- NULL
SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL

head(SeuratTyser_subForMerge_Krup_Tyser_Pijuan@meta.data)

	#Krup
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Batch <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Library.ident
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Dataset <- "Krup.2023"
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Timepoint <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$gem.group
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Old.groups <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$seurat_clusters

v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$gem.group <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$percent.mt <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$Library.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$SCT_snn_res.1 <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data$seurat_clusters <- NULL

head(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan@meta.data)

	# Pijuan
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Batch <- SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$sequencing.batch
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Dataset <- "Pijuan.2019"
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Timepoint <- SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$stage
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Old.groups <- SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$celltype

SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$cell <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$barcode <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$sample <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$pool <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$stage <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$sequencing.batch <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$theiler <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$doub.density <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$doublet <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$cluster <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$cluster.sub <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$cluster.stage <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$cluster.theiler <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$stripped <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$celltype <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$colour <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$sizeFactor <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Limited_Hh_score_1 <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Me4_1 <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Me7_1 <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Me5_1 <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$Me6_1 <- NULL
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL

head(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan@meta.data)

############################
## First quality controls ##
############################
Vln_Krup <- VlnPlot(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Tyser <- VlnPlot(SeuratTyser_subForMerge_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Pijuan <- VlnPlot(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")

plot_Krup1 <- FeatureScatter(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Tyser1 <- FeatureScatter(SeuratTyser_subForMerge_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Pijuan1 <- FeatureScatter(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(file="./1-FirstSeparate-QC_Plots_rawMerge2.pdf", width=9, height=10)
print(Vln_Krup)
print(Vln_Tyser)
print(Vln_Pijuan)
print(plot_Krup1)
print(plot_Tyser1)
print(plot_Pijuan1)
dev.off()

####################
## Merge datasets ##
####################

Merge2_Krup_Tyser_Pijuan <- merge(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan,c(SeuratTyser_subForMerge_Krup_Tyser_Pijuan,SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan))

VarFeatures <- rbind(data.frame(x=VariableFeatures(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan)),data.frame(x=VariableFeatures(SeuratTyser_subForMerge_Krup_Tyser_Pijuan)),data.frame(x=VariableFeatures(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan)))
VarFeatures_List <- unique(VarFeatures$x)

VariableFeatures(Merge2_Krup_Tyser_Pijuan) <- VarFeatures_List

# 39,711 cells - 26,049 SCT features / 6,086 variable features

saveRDS(Merge2_Krup_Tyser_Pijuan, file = "./Merge2_Krup_Tyser_Pijuan.RDS")
#Merge2_Krup_Tyser_Pijuan <- readRDS(file = "./Merge2_Krup_Tyser_Pijuan.RDS")

############################
## Merge quality controls ##
############################
Vln_Merge <- VlnPlot(Merge2_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, group.by = "Dataset")
Vln_MergeSCT <- VlnPlot(Merge2_Krup_Tyser_Pijuan, features = c("nFeature_SCT", "nCount_SCT"), ncol = 3, group.by = "Dataset")
plot_Merge <- FeatureScatter(Merge2_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_MergeSCT <- FeatureScatter(Merge2_Krup_Tyser_Pijuan, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")

pdf(file="./2-Merge2/1-QC_Plots_Merge2.pdf", width=9, height=10)
print(Vln_Merge)
print(Vln_MergeSCT)
print(plot_Merge)
print(plot_MergeSCT)
dev.off()

###################################################
## PCA, Harmony Integration, UMAP and clustering ##
###################################################

Merge2_Krup_Tyser_Pijuan <- RunPCA(Merge2_Krup_Tyser_Pijuan, verbose = TRUE)

pdf(file="./2-Merge2/Merge2_Krup_Tyser_Pijuan_ElbowPlot-PCA.pdf", width=9, height=10)
ElbowPlot(Merge2_Krup_Tyser_Pijuan, ndims = 50, reduction = "pca")
dev.off()
# 30 components

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
Merge2_Krup_Tyser_Pijuan_harmony2 <- RunHarmony(Merge2_Krup_Tyser_Pijuan, 
				group.by.vars = c("Batch","Dataset"), 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

print("Integration / Batch effect correction done!")

# UMAP calculation and clustering
Merge2_Krup_Tyser_Pijuan_harmony2 <- RunUMAP(Merge2_Krup_Tyser_Pijuan_harmony2, dims = 1:30, reduction = "harmony", verbose = TRUE)
Merge2_Krup_Tyser_Pijuan_harmony2 <- FindNeighbors(Merge2_Krup_Tyser_Pijuan_harmony2,reduction = "harmony", dims = 1:30, verbose = TRUE)
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- FindClusters(Merge2_Krup_Tyser_Pijuan_harmony2, verbose = TRUE, resolution = 1.0)

pdf(file="./2-Merge2/UMAP_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T, raster = FALSE))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Dataset", raster = FALSE))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Batch", raster = FALSE))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "orig.ident", raster = FALSE)+ NoLegend())
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Timepoint", raster = FALSE))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE)+ NoLegend())
dev.off()


# Add Log count per million normalization (RNA slot) for gene expression visualization and analysis
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- NormalizeData(
  Merge2_Krup_Tyser_Pijuan_harmony2_res1.0,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1000000,
  verbose = TRUE
)

saveRDS(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
#Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

# Homogeneate timestamps
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge2_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(TimeStamp = case_when(
    Timepoint %in% "E60" ~ "E6.0",
	Timepoint %in% "lateE60" ~ "lateE6.0",
	Timepoint %in% "E65" ~ "E6.5",
	Timepoint %in% "lateE65" ~ "lateE6.5",
	Timepoint %in% "E70" ~ "E7.0",
	Timepoint %in% "earlyE75" ~ "earlyE7.5",
	Timepoint %in% "lateE75" ~ "lateE7.5",
	Timepoint %in% "E775" ~ "E7.75",
	TRUE ~ Timepoint
  ))

# Time Stamp UMAP
pdf(file="./2-Merge2/TimeStamp-UMAP_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
for(Day in unique(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimeStamp)){
cells <- WhichCells(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimeStamp ==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", group.by = "TimeStamp", cells.highlight=Cells, cols.highlight=c("#4472C4"),sizes.highlight = 0.6, pt.size = 0.4, cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

print("Time Stamp UMAP done!")

#######################################################################
## Gene expression analysis - need to switch from SCT to RNA for CPM ##
#######################################################################
TmpObj_GeneExpr <- Merge2_Krup_Tyser_Pijuan_harmony2_res1.0
DefaultAssay(TmpObj_GeneExpr) <- "RNA"

# Dot plot and Feature plots QCs
gene_list <- c("Fgf5","Isl1","Osr1","Tbx1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Asb2","Nkx2-5","Tbx5","Gata4","Gata5","Gata6","Gli1","Foxf1","Hhip","Ptch1","Fgf8","Fgf10","Mab21l2","Hand1","Smarcd3","Ai14-transgene","eGFP-transgene","Myl2","Hcn4","Hopx")
pdf(file = "./2-Merge2/MarkerGenes.v1_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(TmpObj_GeneExpr, features = i) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

saveRDS(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
#Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

################################
## Top markers identification ##
################################
# Find markers
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- PrepSCTFindMarkers(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0)
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.markers <- FindAllMarkers(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6,recorrect_umi = FALSE)
write.csv(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.markers, file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0_all-markers.csv")

# Marker heatmap
Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = "./2-Merge2/TopMarkers_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 11, height = 11)
print(DoHeatmap(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, features = top15$gene) + NoLegend())
dev.off()

pdf(file = "./2-Merge2/DotPlot_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 15, height = 40)
print(DotPlot(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

###################################################
## SingleR - Mouse Gastrulation Atlas annotation ##
###################################################
library(SingleR)

Merge2_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./2-Merge2/Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
Merge2_Krup_Tyser_Pijuan_harmony2_resTest <- FindClusters(Merge2_Krup_Tyser_Pijuan_harmony2_res1.0, verbose = TRUE, resolution = 4)

Merge2_Krup_Tyser_Pijuan_harmony2_resTest <- PrepSCTFindMarkers(Merge2_Krup_Tyser_Pijuan_harmony2_resTest)

ref <- readRDS("./4-singlecell_reference_datasets/RefEmbryoAtlasData.rds")
ref <- ref[, ref$doublet == FALSE] # remove doublets
ref <- ref[, !is.na(ref$celltype)] # remove empty cell types
ref_stage <- ref[, ref$stage %in% c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")] # keep only stage e8.25 - Keep the time stage you want
ref_stage <- scuttle::logNormCounts(ref_stage) # singleR needs log normalized reference
sce_Merge2_Krup_Tyser_Pijuan <- as.SingleCellExperiment(Merge2_Krup_Tyser_Pijuan_harmony2_resTest)

# Convert symbol to ensembl
require(org.Mm.eg.db)
ens <- mapIds(org.Mm.eg.db,
	        keys = rownames(sce_Merge2_Krup_Tyser_Pijuan),
		  column = 'ENSEMBL',
		  keytype = 'SYMBOL')
all(rownames(sce_Merge2_Krup_Tyser_Pijuan) == names(ens))

keep <- !is.na(ens)
ens <- ens[keep]
sce_Merge2_Krup_Tyser_Pijuan <- sce_Merge2_Krup_Tyser_Pijuan[keep,]
rownames(sce_Merge2_Krup_Tyser_Pijuan) <- ens

# Run single R
mouse_gastrulation_label <- SingleR(test=sce_Merge2_Krup_Tyser_Pijuan, ref = ref_stage, labels = ref_stage$celltype, clusters = sce_Merge2_Krup_Tyser_Pijuan$seurat_clusters)
saveRDS(mouse_gastrulation_label, file ="./2-Merge2/SingleR_MouseGastrulationAtlas/Res4/mouse_gastrulation_label.RDS")

Merge2_Krup_Tyser_Pijuan_harmony2_resTest$mouse_gastrulation_label <- factor(Merge2_Krup_Tyser_Pijuan_harmony2_resTest$seurat_clusters, labels = mouse_gastrulation_label$labels)
saveRDS(Merge2_Krup_Tyser_Pijuan_harmony2_resTest, file = "./2-Merge2/SingleR_MouseGastrulationAtlas/Res4/Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR.RDS")
#Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR <- readRDS(file = "./2-Merge2/SingleR_MouseGastrulationAtlas/Res4/Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR.RDS")

pdf(file="./2-Merge2/SingleR_MouseGastrulationAtlas/Res4/UMAP_Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR.pdf", width=9, height=10)
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_resTest, reduction = "umap", group.by = "mouse_gastrulation_label"))
print(DimPlot(Merge2_Krup_Tyser_Pijuan_harmony2_resTest, reduction = "umap", group.by = "mouse_gastrulation_label",label = T)+NoLegend())
dev.off()

pdf(file="./2-Merge2/SingleR_MouseGastrulationAtlas/Res4/UMAP_Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR-plotScoreHeatmap.pdf", width=9, height=10)
print(plotScoreHeatmap(mouse_gastrulation_label, clusters=rownames(mouse_gastrulation_label)))
dev.off()

#===============================================================================================================================================#
######################################################========      Subset 2      ========#######################################################
#===============================================================================================================================================#
# SUBSET 2 preparation - need to remove clusters: 20,21,23,28,31,32,33,38,43 (from Merge2_Krup_Tyser_Pijuan_harmony2_res4.0)
# Need to obtain cell names to remove them from individual dataset followed by separate SCT before merging/Harmony integration
CellsToRemove2 <- as.data.frame(Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR@meta.data[Merge2_Krup_Tyser_Pijuan_harmony2_res4_SlingleR@meta.data$seurat_clusters %in% c(20,21,23,28,31,32,33,38,43),c("seurat_clusters","Dataset","mouse_gastrulation_label")])

CellsToRemove2_Krup <- rownames(CellsToRemove2[CellsToRemove2$Dataset == "Krup.2023",])
write.csv(CellsToRemove2_Krup, file = "./2-Merge2/PrepSub2_CellsToRemove_Krup.csv")

CellsToRemove2_Tyser <- rownames(CellsToRemove2[CellsToRemove2$Dataset == "Tyser.2021",])
write.csv(CellsToRemove2_Tyser, file = "./2-Merge2/PrepSub2_CellsToRemove_Tyser.csv")

CellsToRemove2_Pijuan <- rownames(CellsToRemove2[CellsToRemove2$Dataset == "Pijuan.2019",])
write.csv(CellsToRemove2_Pijuan, file = "./2-Merge2/PrepSub2_CellsToRemove_Pijuan.csv")

#===============================================================================================================================================#
######################################################========      Merge 3      ========#######################################################
#===============================================================================================================================================#
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./Krup_GSE208153_scRNAseq/0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan.RDS")
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./3-PublicDatasets/8-TyserJCF_2021/Integration/SeuratTyser_subForMerge3_Krup_Tyser_Pijuan.RDS")
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./MouseGastrulationData/SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan.RDS")

######################################
## Data cleaning and homogenization ##
######################################
	# Tyser
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$Batch <- SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$batch
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$Dataset <- "Tyser.2021"
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$Timepoint <- SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$stage
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$Old.groups <- SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$clusterAnn

SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$batch <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$sampleID <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$date <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$plate <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$well <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$stage <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$side <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$ENA <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$cluster <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$clusterAnn <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$clusterCol <- NULL
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL

head(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan@meta.data)

	#Krup
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Batch <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Library.ident
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Dataset <- "Krup.2023"
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Timepoint <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$gem.group
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Old.groups <- v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$seurat_clusters

v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$gem.group <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$percent.mt <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$Library.ident <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$SCT_snn_res.1 <- NULL
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data$seurat_clusters <- NULL

head(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan@meta.data)

	# Pijuan
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Batch <- SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$sequencing.batch
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Dataset <- "Pijuan.2019"
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Timepoint <- SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$stage
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Old.groups <- SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$celltype

SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$cell <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$barcode <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$sample <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$pool <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$stage <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$sequencing.batch <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$theiler <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$doub.density <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$doublet <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$cluster <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$cluster.sub <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$cluster.stage <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$cluster.theiler <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$stripped <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$celltype <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$colour <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$sizeFactor <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Limited_Hh_score_1 <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Me4_1 <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Me7_1 <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Me5_1 <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$Me6_1 <- NULL
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data$old.ident <- NULL

head(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan@meta.data)

############################
## First quality controls ##
############################
Vln_Krup <- VlnPlot(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Tyser <- VlnPlot(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")
Vln_Pijuan <- VlnPlot(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident")

plot_Krup1 <- FeatureScatter(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Tyser1 <- FeatureScatter(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_Pijuan1 <- FeatureScatter(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(file="./1-FirstSeparate-QC_Plots_rawMerge3.pdf", width=9, height=10)
print(Vln_Krup)
print(Vln_Tyser)
print(Vln_Pijuan)
print(plot_Krup1)
print(plot_Tyser1)
print(plot_Pijuan1)
dev.off()

####################
## Merge datasets ##
####################
Merge3_Krup_Tyser_Pijuan <- merge(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan,c(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan,SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan))

VarFeatures <- rbind(data.frame(x=VariableFeatures(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan)),data.frame(x=VariableFeatures(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan)),data.frame(x=VariableFeatures(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan)))
VarFeatures_List <- unique(VarFeatures$x)

VariableFeatures(Merge3_Krup_Tyser_Pijuan) <- VarFeatures_List

# 35,617 cells - 25,662 SCT features / 6,212 variable features

saveRDS(Merge3_Krup_Tyser_Pijuan, file = "./Merge3_Krup_Tyser_Pijuan.RDS")
#Merge3_Krup_Tyser_Pijuan <- readRDS(file = "./Merge3_Krup_Tyser_Pijuan.RDS")


############################
## Merge quality controls ##
############################
Vln_Merge <- VlnPlot(Merge3_Krup_Tyser_Pijuan, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, group.by = "Dataset")
Vln_MergeSCT <- VlnPlot(Merge3_Krup_Tyser_Pijuan, features = c("nFeature_SCT", "nCount_SCT"), ncol = 3, group.by = "Dataset")
plot_Merge <- FeatureScatter(Merge3_Krup_Tyser_Pijuan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_MergeSCT <- FeatureScatter(Merge3_Krup_Tyser_Pijuan, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")

pdf(file="./3-Merge3/1-QC_Plots_Merge3.pdf", width=9, height=10)
print(Vln_Merge)
print(Vln_MergeSCT)
print(plot_Merge)
print(plot_MergeSCT)
dev.off()

###################################################
## PCA, Harmony Integration, UMAP and clustering ##
###################################################
Merge3_Krup_Tyser_Pijuan <- RunPCA(Merge3_Krup_Tyser_Pijuan, verbose = TRUE)

pdf(file="./3-Merge3/Merge3_Krup_Tyser_Pijuan_ElbowPlot-PCA.pdf", width=9, height=10)
ElbowPlot(Merge3_Krup_Tyser_Pijuan, ndims = 50, reduction = "pca")
dev.off()
# 30 components

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
Merge3_Krup_Tyser_Pijuan_harmony2 <- RunHarmony(Merge3_Krup_Tyser_Pijuan, 
				group.by.vars = c("Batch","Dataset"), 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

print("Integration / Batch effect correction done!")

# UMAP calculation and clustering
Merge3_Krup_Tyser_Pijuan_harmony2 <- RunUMAP(Merge3_Krup_Tyser_Pijuan_harmony2, dims = 1:30, reduction = "harmony", verbose = TRUE)
Merge3_Krup_Tyser_Pijuan_harmony2 <- FindNeighbors(Merge3_Krup_Tyser_Pijuan_harmony2,reduction = "harmony", dims = 1:30, verbose = TRUE)
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- FindClusters(Merge3_Krup_Tyser_Pijuan_harmony2, verbose = TRUE, resolution = 1.0)

pdf(file="./3-Merge3/UMAP_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T, raster = FALSE))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Dataset", raster = FALSE))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Batch", raster = FALSE))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "orig.ident", raster = FALSE)+ NoLegend())
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Timepoint", raster = FALSE))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Old.groups", raster = FALSE)+ NoLegend())
dev.off()


# Add Log count per million normalization (RNA slot) for gene expression visualization and analysis
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- NormalizeData(
  Merge3_Krup_Tyser_Pijuan_harmony2_res1.0,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1000000,
  verbose = TRUE
)

saveRDS(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
#Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

# Homogeneate timestamps
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(TimeStampMerged = case_when(
    TimeStamp %in% "lateE6.0" ~ "E6.5",
	TimeStamp %in% "lateE6.5" ~ "E6.75",
	TimeStamp %in% "earlyE7.5" ~ "E7.25",
	TimeStamp %in% "lateE7.5" ~ "E7.75",
	TimeStamp %in% "-1" ~ "E7.25",
	TimeStamp %in% "0" ~ "E7.5",
	TimeStamp %in% "1" ~ "E7.75",
	TimeStamp %in% "2" ~ "E8.0",
	TimeStamp %in% "3" ~ "E8.25",
	TimeStamp %in% "LHT" ~ "E8.5",
	TRUE ~ TimeStamp
  ))

# Time Stamp UMAP
pdf(file="./3-Merge3/TimeStampMerged-UMAP_Merge2_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width=9, height=10)
for(Day in unique(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimeStampMerged)){
cells <- WhichCells(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimeStampMerged ==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", group.by = "TimeStampMerged", cells.highlight=Cells, cols.highlight=c("#4472C4"),sizes.highlight = 0.6, pt.size = 0.4, cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

print("Time Stamp UMAP done!")

#######################################################################
## Gene expression analysis - need to switch from SCT to RNA for CPM ##
#######################################################################
TmpObj_GeneExpr <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
DefaultAssay(TmpObj_GeneExpr) <- "RNA"

# Dot plot and Feature plots QCs
gene_list <- c("Fgf5","Isl1","Osr1","Tbx1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Asb2","Nkx2-5","Tbx5","Gata4","Gata5","Gata6","Gli1","Foxf1","Hhip","Ptch1","Fgf8","Fgf10","Nrp1","Mab21l2","Hand1","Hand2","Hey2","Smarcd3","Ai14-transgene","eGFP-transgene","Myl2","Hcn4","Hopx","Smo")
pdf(file = "./3-Merge3/MarkerGenes.v1_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(TmpObj_GeneExpr, features = i) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

saveRDS(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
#Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

################################
## Top markers identification ##
################################
# Find markers
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- PrepSCTFindMarkers(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.markers <- FindAllMarkers(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6,recorrect_umi = FALSE)
write.csv(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.markers, file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_all-markers.csv")

# Marker heatmap
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = "./3-Merge3/TopMarkers_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 11, height = 11)
print(DoHeatmap(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, features = top15$gene) + NoLegend())
dev.off()

pdf(file = "./3-Merge3/DotPlot_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.markers.pdf", width = 15, height = 40)
print(DotPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

######################################
## UCell function Hh score # R/4.3.1 #
######################################
library(Seurat)
library(UCell)
library(ggplot2)
library(circlize)
library(ggh4x)

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- PrepSCTFindMarkers(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)

markers <- list(Limited_Hh_list = c("Ptch1","Hhip","Gli1"))

SeuratObj <- AddModuleScore_UCell(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, markers, assay = "SCT")
SeuratObj <- SmoothKNN(SeuratObj, signature.names = paste0(names(markers),"_UCell"), reduction="harmony", k=10)

UCell_score <- data.frame(UCell = SeuratObj@meta.data$Limited_FHF_list4.2_UCell_kNN, row.names = rownames(SeuratObj@meta.data))
UCell_score2 <- as.data.frame(SeuratObj@meta.data[,paste0(names(markers),"_UCell")])
write.csv(UCell_score, file = "./3-Merge3/UCell/Limited_Hh_list_UCell_kNN_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

