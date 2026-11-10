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
print(DimPlot(Merge_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T,  group.by = "Dataset", cols = c("Krup.2023"="#1E88E5","Pijuan.2019"="#FFC107","Tyser.2021"="#D81B60"), raster = FALSE))
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

## Mixed mesoderm-specific analysis
	# Find global markers
dataset <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
Idents(dataset) <- dataset@meta.data$Annotation
dataset.markers <- FindAllMarkers(dataset, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5, recorrect_umi = FALSE)
write.csv(dataset.markers, file = "./MixedMesoderm_Markers/4-AnnotationCluster_all-markers_min0.5logfc0.5.csv")

	# Global marker heatmap
dataset.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = "./MixedMesoderm_Markers/4-AnnotationCluster_all-markers_min0.5logfc0.5.pdf", width = 15, height = 40)
print(DotPlot(dataset, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 5) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

	# Find pairwise markers
markers <- FindMarkers(
  object = dataset,
  ident.1 = "Mixed mesoderm",
  ident.2 = "Mesenchyme", # test with all individual cell type
  assay = "SCT",          
  slot = "data",          
  test.use = "wilcox",    
  logfc.threshold = 0.6,  
  min.pct = 0.5,         
  min.diff.pct = 0.2,     
  min.cells.group = 10,   
  only.pos = TRUE,        # focus on positive markers of "Mixed mesoderm"
  max.cells.per.ident = 2000
)
write.csv(markers, file = "./MixedMesoderm_Markers/FindMarkers_min0.5logfc0.6_MixedMesoderm-vs-Mesenchyme.csv")


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

UCell_score <- data.frame(UCell = SeuratObj@meta.data$Limited_Hh_list_UCell_kNN, row.names = rownames(SeuratObj@meta.data))
write.csv(UCell_score, file = "./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

##########################################################################
## Early/Late mixed mesoderm analysis and comparison with Lescroart.2014 #
##########################################################################
dataset <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0

## Early/Late mixed mesoderm analysis

Subset_MixedMonly <- subset(x = dataset, subset = Annotation %in% "Mixed mesoderm", invert=FALSE)
Subset_MixedMonly_SubTime <- subset(x = Subset_MixedMonly, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.25","E7.5"), invert=FALSE)

library(DESeq2)
set.seed(123)

	# Assign cells randomly to one of 3 groups
Subset_MixedMonly_SubTime$pseudoGroup <- sample(1:3, nrow(Subset_MixedMonly_SubTime@meta.data), replace = T)

Subset_MixedMonly_SubTime@meta.data <- Subset_MixedMonly_SubTime@meta.data %>%
  mutate(Condition = case_when(
	TimeStampMerged %in% "E6.5" ~ "Early",
	TimeStampMerged %in% "E6.75" ~ "Early",
	TimeStampMerged %in% "E7.25" ~ "Late",
	TimeStampMerged %in% "E7.5" ~ "Late",
  ))
saveRDS(Subset_MixedMonly_SubTime, file = "./Pseudobulk_comparisons/5g-Subset_MixedMonly_SubTime.RDS")

	# Aggregate expression - make sure to always include pseudoGroup (or whatever group) in the group by.
pseudo_Subset_MixedMonly_SubTime_counts <- AggregateExpression(Subset_MixedMonly_SubTime, assay="RNA", slot="counts",return.seurat = F, group.by = c("Condition", "pseudoGroup"))
pseudo_Subset_MixedMonly_SubTime_counts <- as.matrix(pseudo_Subset_MixedMonly_SubTime_counts$RNA)

	# Remove mitochondrial and ribosomal genes
keep <- grep("^mt-|^Rpl|^Rps|^Mrpl*", rownames(pseudo_Subset_MixedMonly_SubTime_counts), invert = T)
pseudo_Subset_MixedMonly_SubTime_counts_mt_rb_removed <- pseudo_Subset_MixedMonly_SubTime_counts[keep, ]

dim(pseudo_Subset_MixedMonly_SubTime_counts) ; dim(pseudo_Subset_MixedMonly_SubTime_counts_mt_rb_removed)

saveRDS(pseudo_Subset_MixedMonly_SubTime_counts, file = "./Pseudobulk_comparisons/5g-pseudo_Subset_MixedMonly_SubTime_counts.RDS")

	# Create meta data for deseq
metad <- data.frame(samples = colnames(pseudo_Subset_MixedMonly_SubTime_counts)) %>%
  mutate(group = str_split_fixed(colnames(pseudo_Subset_MixedMonly_SubTime_counts), "_",2)[,2] %>% as.factor(),
		 condition = str_split_fixed(colnames(pseudo_Subset_MixedMonly_SubTime_counts), "_",2)[,1] %>% as.factor())
 rownames(metad) <- metad$samples

metad$condition <- relevel(metad$condition, ref="Early")

	# remove lowely expressed genes with edgeR
library(edgeR)

mm <- model.matrix(~ group + condition, data=metad)
keep <- filterByExpr(pseudo_Subset_MixedMonly_SubTime_counts_mt_rb_removed, design = mm)
pseudo_Subset_MixedMonly_SubTime_counts_filtered <- pseudo_Subset_MixedMonly_SubTime_counts_mt_rb_removed[keep, ]
dim(pseudo_Subset_MixedMonly_SubTime_counts_mt_rb_removed)
dim(pseudo_Subset_MixedMonly_SubTime_counts_filtered)

dds <- DESeqDataSetFromMatrix(pseudo_Subset_MixedMonly_SubTime_counts_filtered, metad, design = ~ group + condition)
dds <- DESeq(dds)

	#saves
saveRDS(dds, file = "./Pseudobulk_comparisons/5g-MixedMesoderm_EarlyLate_dds.RDS")

	# PCA
library(FactoMineR)
library(factoextra)

vsd <- vst(dds)
pca_result <- PCA(t(assay(vsd)), graph = FALSE)
pca_data <- data.frame(samples = rownames(pca_result$ind$coord), pca_result$ind$coord)
pca_data <- left_join(pca_data, metad, by = "samples")

	# Plotting the PCA results
pdf(file="./Pseudobulk_comparisons/5g-MixedMesoderm_EarlyLate_PCA.pdf", width=9, height=10)
pc1_var <- round(pca_result$eig[1, 2], 1)
pc2_var <- round(pca_result$eig[2, 2], 1)
ggplot(pca_data, aes(x = Dim.1, y = Dim.2, color = condition, shape = group)) +
  geom_point(size = 3) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  )
dev.off()
  
	# DESeq2 results 
levels(dds$condition)
contrasts <- list(MixedMeso = c("Late", "Early"))

res_deseq<- list()
for (celltype in names(contrasts)) {
  contrast = contrasts[[celltype]]
  res = results(dds, contrast = c("condition",contrast[1], contrast[2]), alpha=0.05)
  res_deseq[[celltype]] <- res
}
write.csv(res_deseq, file = "./Pseudobulk_comparisons/5g-MixedMesoderm_EarlyLate_DESeq_results.csv")

	# Volcano plot
source("./gg_volcano4.R")

pdf(file="./Pseudobulk_comparisons/5g-MixedMesoderm_EarlyLate_Volcano.pdf", width=5, height=5)
gg_volcano4(as.data.frame(res_deseq$MixedMeso),
			abs_lfc= 1 , 
			p.cutoff = 0.05,
			up.color="#C2202C",
			down.color="#0070C0",
			unSig.color='grey88', 
			dotSize=0.5, 
			alpha=1,    
			theme_style = "minimal",
			scale_axes = TRUE, # TRUE for auto-scaling, FALSE for manual
			xlim = NULL,  
			ylim = NULL, 
			count_position = "bottom",
			use_asinh = TRUE, # TRUE to use asinh transformation on y-axis
			n_labels_up = 10,
			n_labels_down = 10,
			label_TF_only = FALSE,
			align_labels = FALSE
			)
gg_volcano4(as.data.frame(res_deseq$MixedMeso),
			abs_lfc= 1 , 
			p.cutoff = 0.05,
			up.color="#C2202C",
			down.color="#0070C0",
			unSig.color='grey88', 
			dotSize=0.5, 
			alpha=1,    
			theme_style = "minimal",
			scale_axes = TRUE, # TRUE for auto-scaling, FALSE for manual
			xlim = NULL,
			ylim = NULL,
			count_position = "bottom",
			use_asinh = TRUE, # TRUE to use asinh transformation on y-axis
			n_labels_up = 10,
			n_labels_down = 10,
			label_TF_only = TRUE,
			align_labels = TRUE
			)
dev.off()

	# GO term analysis
library(clusterProfiler)
library(org.Mm.eg.db)

	# Load gene annotation
eg <- read_csv("./Mus_musculus.gene_info_edit.csv")
eg$GeneID <- as.character(eg$GeneID)

run_go_local <- function(test.vector, univr_vector, ont = "BP") {
  test.ids  <- subset(eg, Symbol %in% test.vector  | splt_syn %in% test.vector)$GeneID  %>% unique()
  univr.ids <- subset(eg, Symbol %in% univr_vector | splt_syn %in% univr_vector)$GeneID %>% unique()

  goPATH <- enrichGO(
    gene          = test.ids,
    OrgDb         = "org.Mm.eg.db",
    universe      = univr.ids,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05
  )
  goPATH <- setReadable(goPATH, OrgDb = "org.Mm.eg.db")
  return(goPATH)
}

	# Extract gene lists from DESeq2 results
res_df <- read.csv("./Pseudobulk_comparisons/5g-MixedMesoderm_EarlyLate_DESeq_results.csv") %>%
  rename(gene_id = X)

universe_genes <- res_df$gene_id
up_genes       <- subset(res_df, MixedMeso.padj < 0.05 & MixedMeso.log2FoldChange >  1)$gene_id
down_genes     <- subset(res_df, MixedMeso.padj < 0.05 & MixedMeso.log2FoldChange < -1)$gene_id

stopifnot(length(up_genes) > 5, length(down_genes) > 5)

	# Run GO
go_up   <- run_go_local(up_genes,   universe_genes)
go_down <- run_go_local(down_genes, universe_genes)

	# Plots
out_dir <- "./Pseudobulk_comparisons/"

plot_go <- function(go_obj, title, color) {
  dotplot(go_obj, showCategory = 20) +
    ggtitle(title) +
    scale_color_gradient(low = color, high = "grey80") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3, linetype = "dashed"),
      panel.grid.minor   = element_blank(),
      legend.text        = element_text(size = 8),
      legend.title       = element_text(size = 8),
      axis.line.x        = element_line(color = "black"),
      axis.line.y        = element_line(color = "black"),
      axis.ticks         = element_line(color = "black"),
      plot.margin        = margin(5, 5, 5, 80)
    ) +
    scale_size(range = c(3, 10)) +
    scale_x_continuous(expand = expansion(mult = c(0.5, 0.15)))
}

pdf(file.path(out_dir, "5g-MixedMesoderm_EarlyLate_GO.pdf"), width = 9, height = 8)
print(plot_go(go_up,   "Upregulated GO BP\nLate vs Early",   "#C2202C"))
print(plot_go(go_down, "Downregulated GO BP\nLate vs Early", "#0070C0"))
dev.off()

	# Save tables
write.csv(as.data.frame(go_up),   file.path(out_dir, "5g-MixedMesoderm_EarlyLate_GO_up.csv"))
write.csv(as.data.frame(go_down), file.path(out_dir, "5g-MixedMesoderm_EarlyLate_GO_down.csv"))

## Comparison to Lescroart.2014 analysis
library(ggrepel)

	# List abs(FC) > 1.5 + Genes on plot Fig5k-m
E6.5_genes <- c("Aplnr","Mixl1","Frzb","Meis2","Fn1","Mesp1","Dll3","Cdh11","Gata6","Dlc1",
           "Fgf3","Sp5","Pcdh7","Adam19","Hoxb1","Gas1","Lhx1","Peg3","Gna14","Pdgfra",
           "Snai1","Wnt5a","Epha4","T","Ifitm1","Tbx3","Klhl6","Bmp7","Cyp26a1","Nin",
           "Cxcr4","Flk1","Nfatc1","Gata4","Prkd1","Rab8b","Rnd3","Cdh2","Ptprj","Zfx",
           "Spred1","Fgf10","Ccnd2","Vangl1","Mesp2","Tbx6","Pdlim5","Otx1","Aph1b","Smad1",
           "L1cam","Chst2","Man1c1","Evx1","Mmp14","Irx5","Tgfb1","Cdh4","Sall3","Prtg",
           "Usp3","Asb4","Arf1","Wwp1","Cdx1","Cdc42ep4","Sema5a","Birc6","Whsc1l1","Vcan",
           "Odz4","Wnt3","Lefty2","Cbln1","Anpep","Rimbp2","Efna3","Vldlr","Slit3","Id1",
           "Rftn1","Cachd1","Atp11c","Hand1","Dlk1","Pafah1b1","Zfp423","Msgn1","Apln","Efnb3",
           "Hs3st3b1","Six2","Rreb1","Wnt2b","Dact1","Fstl1","Rasgrp3","Gpr50","Braf","Adcyap1r1",
           "Tiparp","Dock11","Itga5","Agtrap","Nrp2","Dclk1","Phlda2","Ppic","Psmd9","Wls",
           "Fzd1","Axin2","Atl1","Itgav","Chst7","Man1a2","Bmper","S1pr5","Arl4d","Pcdh18",
           "Nfat5","Zeb2","Glis1","Wars2","Gabpa","Lef1","Phc2","Cask","Chek1","Pbx1",
           "Irx3","Fat3","Atg5",
		   "Eomes","Fgf8","Otx2","Nanog")

E7.5_genes <- c("Zfhx3","Nasp","Cer1","Ptpra","Rgnef","Runx2","Gfra2","Actc1","Rxra","Tdgf1",
           "Egfr","Txnrd2","Gys1","Klf12","Dab2","Meox1","Tnnc1",
		   "Meis1","Hoxb3","Hoxa1","Hoxa3","Foxc2","Aldh1a2","Cited1","Myl7","Hoxb2","Foxc1","Smarcd3")
		   
	# Build plotting data frame
plot_df <- bind_rows(
  res_df %>% filter(gene_id %in% E6.5_genes) %>%
    mutate(group = "E6.5 Mesp1+ (Lescroart 2014)"),
  res_df %>% filter(gene_id %in% E7.5_genes) %>%
    mutate(group = "E7.5 Mesp1+ (Lescroart 2014)")
) %>%
  dplyr::select(gene_id, log2FoldChange = MixedMeso.log2FoldChange, group)

# Plot
pos <- position_jitter(width = 0.1, seed = 42)
counts_df <- plot_df %>%
  group_by(group) %>%
  summarise(
    n_up   = sum(log2FoldChange >  1, na.rm = TRUE),
    n_down = sum(log2FoldChange < -1, na.rm = TRUE),
    .groups = "drop"
  )

pdf(file.path(out_dir, "5g-MixedMesoderm_Lescroart_violin_absFCdiff_sup1.5_PlusGenePlot.pdf"), width = 5, height = 5)
ggplot(plot_df, aes(x = group, y = log2FoldChange, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_point(aes(color = group), position = pos, size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept =  1, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_text(data = counts_df, aes(x = group, label = paste0("log2FC > 1 = ", n_up),   color = group),
            y = max(plot_df$log2FoldChange, na.rm = TRUE)+2, vjust = -0.5, size = 3, inherit.aes = FALSE) +
  geom_text(data = counts_df, aes(x = group, label = paste0("log2FC < 1 = ", n_down), color = group),
            y = min(plot_df$log2FoldChange, na.rm = TRUE)-1, vjust =  1.5, size = 3, inherit.aes = FALSE) +
  geom_text_repel(data = filter(plot_df, abs(log2FoldChange) > 1),
                  aes(label = gene_id, color = group),
                  position = pos,
                  size = 2.5, max.overlaps = 30,
                  min.segment.length = Inf,
                  segment.size = 0.3) +
  scale_fill_manual(values  = c("E6.5 Mesp1+ (Lescroart 2014)" = "#0070C0", "E7.5 Mesp1+ (Lescroart 2014)" = "#C2202C")) +
  scale_color_manual(values = c("E6.5 Mesp1+ (Lescroart 2014)" = "#0070C0", "E7.5 Mesp1+ (Lescroart 2014)" = "#C2202C")) +
  scale_y_continuous(breaks = sort(unique(c(-1, 0, 1,
    scales::extended_breaks(n = 5)(range(plot_df$log2FoldChange, na.rm = TRUE))))),
    expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = NULL, y = "Log2 Fold Change\n(Late Mixed mesoderm vs Early Mixed mesoderm)") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid     = element_blank(),
        axis.line      = element_line(color = "black"),
        axis.ticks     = element_line(color = "black"))
dev.off()
