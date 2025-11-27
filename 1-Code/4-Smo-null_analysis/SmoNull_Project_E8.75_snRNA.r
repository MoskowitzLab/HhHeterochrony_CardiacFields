#######################################################################################################################################################
################################################                                                       ################################################
################################################ Smo null project - Whole embryo E8.75 - (v11-27-2025) ################################################
################################################                                                       ################################################
#######################################################################################################################################################
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
  
#======================================================================================================================================================#
######################################################========      Data preparation      ========######################################################
#======================================================================================================================================================#

#######################
## Data directories  ##
#######################

# Filtered CellRanger counts
	# Batch1
data_dir1_s1 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/1-Rep1_2samples_1lane_102125/IM-RC-10X-2s-RC12-1-1015/outs/filtered_feature_bc_matrix'
data_dir1_s2 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/1-Rep1_2samples_1lane_102125/IM-RC-10X-2s-RC12-2-1015/outs/filtered_feature_bc_matrix'

	# Batch2
data_dir2_s1 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/2-Rep2_4samples_3lanes_110625/IM-RC-10X-4s-RC13-01/outs/filtered_feature_bc_matrix'
data_dir2_s2 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/2-Rep2_4samples_3lanes_110625/IM-RC-10X-4s-RC13-02/outs/filtered_feature_bc_matrix'
data_dir2_s3 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/2-Rep2_4samples_3lanes_110625/IM-RC-10X-4s-RC13-03/outs/filtered_feature_bc_matrix'
data_dir2_s4 <- '/project/imoskowitz/rcanac/14-SmoNull_snRNA/1-CellRanger/2-Rep2_4samples_3lanes_110625/IM-RC-10X-4s-RC13-04/outs/filtered_feature_bc_matrix'

print(list.files(data_dir1_s1)) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

####################
## Seurat objects ##
####################
Sample1_Control_rep1.data = Read10X(data.dir = data_dir1_s1)
Sample2_SmoNull_rep1.data = Read10X(data.dir = data_dir1_s2)
Sample3_Control_rep2.data = Read10X(data.dir = data_dir2_s1)
Sample4_SmoNull_rep2.data = Read10X(data.dir = data_dir2_s2)
Sample5_Control_rep3.data = Read10X(data.dir = data_dir2_s3)
Sample6_SmoNull_rep3.data = Read10X(data.dir = data_dir2_s4)

Sample1_Control_rep1 = CreateSeuratObject(counts = Sample1_Control_rep1.data)
Sample2_SmoNull_rep1 = CreateSeuratObject(counts = Sample2_SmoNull_rep1.data)
Sample3_Control_rep2 = CreateSeuratObject(counts = Sample3_Control_rep2.data)
Sample4_SmoNull_rep2 = CreateSeuratObject(counts = Sample4_SmoNull_rep2.data)
Sample5_Control_rep3 = CreateSeuratObject(counts = Sample5_Control_rep3.data)
Sample6_SmoNull_rep3 = CreateSeuratObject(counts = Sample6_SmoNull_rep3.data)

# Rename cells
Sample1_Control_rep1 <- RenameCells(object = Sample1_Control_rep1, add.cell.id = "Sample1_Control_rep1")
Sample2_SmoNull_rep1 <- RenameCells(object = Sample2_SmoNull_rep1, add.cell.id = "Sample2_SmoNull_rep1")
Sample3_Control_rep2 <- RenameCells(object = Sample3_Control_rep2, add.cell.id = "Sample3_Control_rep2")
Sample4_SmoNull_rep2 <- RenameCells(object = Sample4_SmoNull_rep2, add.cell.id = "Sample4_SmoNull_rep2")
Sample5_Control_rep3 <- RenameCells(object = Sample5_Control_rep3, add.cell.id = "Sample5_Control_rep3")
Sample6_SmoNull_rep3 <- RenameCells(object = Sample6_SmoNull_rep3, add.cell.id = "Sample6_SmoNull_rep3")

# Add metadata
Sample1_Control_rep1$"Condition" <- "Control"
Sample1_Control_rep1$"Library.ident" <- "Batch1"
Sample1_Control_rep1$"Sample" <- "Sample1"
Sample1_Control_rep1$"Replicate" <- "Rep1"

Sample2_SmoNull_rep1$"Condition" <- "SmoNull"
Sample2_SmoNull_rep1$"Library.ident" <- "Batch1"
Sample2_SmoNull_rep1$"Sample" <- "Sample2"
Sample2_SmoNull_rep1$"Replicate" <- "Rep1"

Sample3_Control_rep2$"Condition" <- "Control"
Sample3_Control_rep2$"Library.ident" <- "Batch2"
Sample3_Control_rep2$"Sample" <- "Sample3"
Sample3_Control_rep2$"Replicate" <- "Rep2"

Sample4_SmoNull_rep2$"Condition" <- "SmoNull"
Sample4_SmoNull_rep2$"Library.ident" <- "Batch2"
Sample4_SmoNull_rep2$"Sample" <- "Sample4"
Sample4_SmoNull_rep2$"Replicate" <- "Rep2"

Sample5_Control_rep3$"Condition" <- "Control"
Sample5_Control_rep3$"Library.ident" <- "Batch2"
Sample5_Control_rep3$"Sample" <- "Sample5"
Sample5_Control_rep3$"Replicate" <- "Rep3"

Sample6_SmoNull_rep3$"Condition" <- "SmoNull"
Sample6_SmoNull_rep3$"Library.ident" <- "Batch2"
Sample6_SmoNull_rep3$"Sample" <- "Sample6"
Sample6_SmoNull_rep3$"Replicate" <- "Rep3"

# Merge samples
E875_SmoNull_Merge <- merge(Sample1_Control_rep1, c(Sample2_SmoNull_rep1,Sample3_Control_rep2,Sample4_SmoNull_rep2,Sample5_Control_rep3,Sample6_SmoNull_rep3))

E875_SmoNull_Merge
table(E875_SmoNull_Merge@meta.data$Sample)

saveRDS(E875_SmoNull_Merge, file = "./E875_SmoNull_Merge_raw.RDS")
#E875_SmoNull_Merge <- readRDS(file = "./E875_SmoNull_Merge_raw.RDS")
#210,779
#Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
#  33974   19276   46350   35500   31459   44220

#======================================================================================================================================================#
######################################################========      Quality controls      ========######################################################
#======================================================================================================================================================#

# Percent mito
E875_SmoNull_Merge[["percent.mt"]] <- PercentageFeatureSet(E875_SmoNull_Merge, pattern = "^mt-")

Idents(E875_SmoNull_Merge) <- E875_SmoNull_Merge@meta.data$Sample
# Visualize QC metrics as merged objects
pdf(file="./1-QCs/QC_Plots.pdf", width=9, height=10)
#print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, split.by = 'Sample'))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0))
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample"))
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample"))
dev.off()

#Test different cutoffs (manually)
E875_SmoNull_MergeSub <- subset(E875_SmoNull_Merge, subset = nCount_RNA > 3000)
E875_SmoNull_MergeSub <- subset(E875_SmoNull_MergeSub, subset = percent.mt < 5)
pdf(file="./1-QCs/QCtest_Plots_Count3000_and_5pMito.pdf", width=9, height=10)
print(VlnPlot(E875_SmoNull_MergeSub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
print(VlnPlot(E875_SmoNull_MergeSub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0))
print(VlnPlot(E875_SmoNull_MergeSub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, split.by = 'Condition'))
print(VlnPlot(E875_SmoNull_MergeSub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0, split.by = 'Condition'))
dev.off()
###

# Subset percent.mt and nCount

# 1. Filter cells by nCount_RNA (UMI counts)
E875_SmoNull_Merge <- subset(E875_SmoNull_Merge, subset = nCount_RNA > 3000)
pdf(file="./1-QCs/QC_Plots_subCount3000.pdf", width=9, height=10)
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Condition"))
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Condition"))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, split.by = 'Condition'))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0, split.by = 'Condition'))
dev.off()

# 2. Filter cells by percent.mito 
E875_SmoNull_Merge <- subset(E875_SmoNull_Merge, subset = percent.mt < 5)
pdf(file="./1-QCs/QC_Plots_subCount3000Mito5.pdf", width=9, height=10)
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Condition"))
print(FeatureScatter(E875_SmoNull_Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Condition"))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, split.by = 'Condition'))
print(VlnPlot(E875_SmoNull_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0, split.by = 'Condition'))
dev.off()

saveRDS(E875_SmoNull_Merge, file = "./E875_SmoNull_Merge_subCount3000Mito5.RDS")
#E875_SmoNull_Merge <- readRDS(file = "./E875_SmoNull_Merge_subCount3000Mito5.RDS")
#168,595

table(E875_SmoNull_Merge@meta.data$Sample)
#Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
#  28470   16833   38791   24110   24973   35418

#========================================================================================================================================================================#
######################################################========      Doublet identification and removal      ========######################################################
#========================================================================================================================================================================#

library(DoubletFinder) #v2.0.3 for compatibility with S4 objects at paramSweep_v3 step

E875_SmoNull_Merge.split <- SplitObject(E875_SmoNull_Merge, split.by = "Sample") 

# Based on the number of called cells by CellRanger: 
# Sample 1 (33,974 cells): 14% doublet
# Sample 2 (19,276 cells): 8% doublet
# Sample 3 (46,350 cells): 19% doublet
# Sample 4 (35,500 cells): 14% doublet
# Sample 5 (31,459 cells): 13% doublet
# Sample 6 (44,220 cells): 18% doublet

PercList <- c(14,8,19,14,13,18)

# 1- Doublet Finder loop
loopLog <- NULL
for(samp in seq(3,length(E875_SmoNull_Merge.split))){
E875_SmoNull_Merge.loop <- NULL

 # Pre-process seurat object with SCT without regression
E875_SmoNull_Merge.loop <- SCTransform(E875_SmoNull_Merge.split[[samp]], verbose = TRUE)
E875_SmoNull_Merge.loop <- RunPCA(E875_SmoNull_Merge.loop, verbose = TRUE)

  # Find significant PCs
  stdv <- E875_SmoNull_Merge.loop[["pca"]]@stdev
  sum.stdv <- sum(E875_SmoNull_Merge.loop[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  E875_SmoNull_Merge.loop <- RunUMAP(E875_SmoNull_Merge.loop, dims = 1:min.pc)
  E875_SmoNull_Merge.loop <- FindNeighbors(object = E875_SmoNull_Merge.loop, dims = 1:min.pc)              
  E875_SmoNull_Merge.loop <- FindClusters(object = E875_SmoNull_Merge.loop, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(E875_SmoNull_Merge.loop, PCs = 1:min.pc, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- E875_SmoNull_Merge.loop@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round((PercList[[samp]]/100) * nrow(E875_SmoNull_Merge.loop@meta.data)) # Change doublet percentage
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  E875_SmoNull_Merge.loop <- doubletFinder_v3(seu = E875_SmoNull_Merge.loop, 
                                   PCs = 1:min.pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj,
								   sct = TRUE)
  metadata <- E875_SmoNull_Merge.loop@meta.data
  colnames(metadata)[13] <- "doublet_finder" # adjust col number
  colnames(metadata)[14] <- "doublet_status" # adjust col number
  E875_SmoNull_Merge.loop@meta.data <- metadata
  
write.csv(E875_SmoNull_Merge.loop@meta.data, file = paste0("./2-DoubletFinder/DoubletData_Sample",samp,"_Doublet-",PercList[[samp]],"perc_E875_SmoNull_Merge.csv"))

# Sample UMAP by doublet
pdf(file=paste0("./2-DoubletFinder/UMAPDoublet_Sample",samp,"_Doublet-",PercList[[samp]],"perc_E875_SmoNull_Merge.pdf"), width=9, height=10)
print(DimPlot(E875_SmoNull_Merge.loop, pt.size = 0.7, reduction = "umap", label = T,  group.by = "doublet_status", raster = FALSE))
dev.off()

# Basic line plot with points
loop <- data.frame(Perc = E875_SmoNull_Merge.loop@meta.data$doublet_finder, Doublet = E875_SmoNull_Merge.loop@meta.data$doublet_status, row.names = rownames(E875_SmoNull_Merge.loop@meta.data))

pdf(file=paste0("./2-DoubletFinder/PlotDoublet_Sample",samp,"_Doublet-",PercList[[samp]],"perc_E875_SmoNull_Merge.pdf"), width=9, height=10)
print(ggplot(data=loop, aes(x=forcats::fct_inorder(rownames(loop)), y=Perc, colour = Doublet)) +
  geom_line()+
  labs(title="Doublet score 14% cutoff",x="Cells", y = "DoubletScore")+
  geom_point()+
  theme(axis.text.x = element_blank()))
dev.off()

loopLog <- c(loopLog, print(paste0("Sample",samp," - ",PercList[[samp]],"% doublet - Done!")))
loopLog <- c(loopLog, print(table(E875_SmoNull_Merge.loop@meta.data$doublet_status)))
}
write.csv(loopLog, file = paste0("./2-DoubletFinder/loopLog2.csv"))

# 2- Flag Doublet cells in the merged object
Samp1 <- read.csv("./2-DoubletFinder/DoubletData_Sample1_Doublet-14perc_E875_SmoNull_Merge.csv")
Samp2 <- read.csv("./2-DoubletFinder/DoubletData_Sample2_Doublet-8perc_E875_SmoNull_Merge.csv")
Samp3 <- read.csv("./2-DoubletFinder/DoubletData_Sample3_Doublet-19perc_E875_SmoNull_Merge.csv")
Samp4 <- read.csv("./2-DoubletFinder/DoubletData_Sample4_Doublet-14perc_E875_SmoNull_Merge.csv")
Samp5 <- read.csv("./2-DoubletFinder/DoubletData_Sample5_Doublet-13perc_E875_SmoNull_Merge.csv")
Samp6 <- read.csv("./2-DoubletFinder/DoubletData_Sample6_Doublet-18perc_E875_SmoNull_Merge.csv")

mergeMeta <- rbind(Samp1,Samp2,Samp3,Samp4,Samp5,Samp6)
rownames(mergeMeta) <- mergeMeta$X
mergeMeta$X <- NULL

E875_SmoNull_Merge@meta.data <- mergeMeta

# Merge
#Doublet Singlet
#  21414  147181

# 3- Remove Doublet cells
E875_SmoNull_Merge_RmDB <- subset(E875_SmoNull_Merge, subset =  doublet_status == "Singlet", invert = FALSE)

#===================================================================================================================================================================#
######################################################========      Normalization and integration      ========######################################################
#===================================================================================================================================================================#
#########################
## Cell cycle scoring  ##
#########################
cc.genes <- read.table("/project/imoskowitz/rcanac/0-R_codes/2-Ressources/Cell-Cycle_genes_scRNA-Seq_MusMusculus.txt", sep="\t", header = T,fill = TRUE)
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[cc.genes$CellCyclePhase == "S",1]
g2m.genes <- cc.genes[cc.genes$CellCyclePhase == "G2/M",1]

# Assign cell cycle scores
E875_SmoNull_Merge_RmDB <- SCTransform(E875_SmoNull_Merge_RmDB, vars.to.regress = c("percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE, conserve.memory = TRUE)
E875_SmoNull_Merge_RmDB <- CellCycleScoring(E875_SmoNull_Merge_RmDB, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = E875_SmoNull_Merge_RmDB@meta.data)

saveRDS(E875_SmoNull_Merge_RmDB, file = "./E875_SmoNull_Merge_RmDB_PostSCTPriorCCcorrect.RDS")
#E875_SmoNull_Merge_RmDB <- readRDS(file = "./E875_SmoNull_Merge_RmDB_PostSCTPriorCCcorrect.RDS")

#########################################
## Normalization, UMAP, Feature plots  ##
#########################################
E875_SmoNull_Merge_RmDB <- SCTransform(E875_SmoNull_Merge_RmDB, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE, conserve.memory = TRUE)
E875_SmoNull_Merge_RmDB_pca <- RunPCA(E875_SmoNull_Merge_RmDB, verbose = TRUE)

pdf(file="./ElbowPlot-PCA_E875_SmoNull_Merge_RmDB_pca.pdf", width=9, height=10)
ElbowPlot(E875_SmoNull_Merge_RmDB_pca, ndims = 50, reduction = "pca")
dev.off()
# 30 components
### Harmony integration by sequencing batch

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
E875_SmoNull_Merge_RmDB_harmonySeqBatch <- RunHarmony(E875_SmoNull_Merge_RmDB_pca, 
				group.by.vars = "Library.ident", 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

# UMAP calculation and clustering
E875_SmoNull_Merge_RmDB_harmonySeqBatch <- RunUMAP(E875_SmoNull_Merge_RmDB_harmonySeqBatch, dims = 1:30, reduction = "harmony", verbose = TRUE)
E875_SmoNull_Merge_RmDB_harmonySeqBatch <- FindNeighbors(E875_SmoNull_Merge_RmDB_harmonySeqBatch,reduction = "harmony", dims = 1:30, verbose = TRUE)
E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0 <- FindClusters(E875_SmoNull_Merge_RmDB_harmonySeqBatch, verbose = TRUE, resolution = 1.0)

# UMAP DimPlots
pdf(file="./UMAP_E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T, raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T,  group.by = "Condition", raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T,  group.by = "Sample", raster = FALSE))
dev.off()

saveRDS(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, file = "./E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0.RDS")
#E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0 <- readRDS(file = "./E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0.RDS")

#Sample1 Sample2 Sample3 Sample4 Sample5 Sample6
#  25144   15650   32772   21233   22312   30070

pdf(file="./UMAP-PerRep_E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T, split.by = "Replicate", group.by = "Condition", raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T, split.by = "Replicate",  group.by = "Sample", raster = FALSE))
dev.off()

############################################################
## SingleR annotation - Extended mouse gastrulation atlas ##
############################################################
library(SingleR)
dataset <- readRDS(file = "./E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0.RDS")

mouse_atlas_extended <- readRDS("./4-singlecell_reference_datasets/EmbryoSCE_extended_atlas.rds")

mouse_atlas_extended  <- mouse_atlas_extended [, mouse_atlas_extended$stage %in% c("E8.5","E8.75","E9.0")]
mouse_atlas_extended <- scuttle::logNormCounts(mouse_atlas_extended) # singleR needs log normalized reference
sce_dataset <- as.SingleCellExperiment(dataset)

# Convert symbol to ensembl
require(org.Mm.eg.db)
ens <- mapIds(org.Mm.eg.db,
	        keys = rownames(sce_dataset),
		  column = 'ENSEMBL',
		  keytype = 'SYMBOL')
all(rownames(sce_dataset) == names(ens))

keep <- !is.na(ens)
ens <- ens[keep]
sce_dataset <- sce_dataset[keep,]
rownames(sce_dataset) <- ens

# Run single R
Extended_mouse_gastrulation_label <- SingleR(test=sce_dataset, ref = mouse_atlas_extended, labels = mouse_atlas_extended$celltype_extended_atlas, clusters = sce_dataset$seurat_clusters)

saveRDS(Extended_mouse_gastrulation_label, file ="./3-SingleR/Extended_mouse_gastrulation_label.RDS")
#Extended_mouse_gastrulation_label <-readRDS("./3-SingleR/Extended_mouse_gastrulation_label.RDS")

dataset$Extended_mouse_gastrulation_label <- factor(dataset$seurat_clusters, labels = Extended_mouse_gastrulation_label$labels)

pdf(file="./3-SingleR/UMAP_E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0_SlingleR.pdf", width=9, height=10)
print(DimPlot(dataset, reduction = "umap", group.by = "Extended_mouse_gastrulation_label",raster=FALSE))
print(DimPlot(dataset, reduction = "umap", group.by = "Extended_mouse_gastrulation_label",label = T,raster=FALSE)+NoLegend())
dev.off()
pdf(file="./3-SingleR/Heatmap_E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0_SlingleR-plotScoreHeatmap.pdf", width=9, height=10)
print(plotScoreHeatmap(Extended_mouse_gastrulation_label, clusters=rownames(Extended_mouse_gastrulation_label)))
dev.off()

saveRDS(dataset, file = "./3-SingleR/E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0_SlingleR.RDS")

#================================================================================================================================================================#
######################################################========      Cardiac lineage refinement      ========######################################################
#================================================================================================================================================================#
## keep clusters from res1.0 -> 2,17,18,20,23,27,30,35 - Cardiopharyngeal progenitors, Mesenchyme, Epicardium, Cardiomyocytes

E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1 <- subset(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0, idents = c(2,17,18,20,23,27,30,35), invert=FALSE)
# 25,618 cells

#########################################
## Normalization, UMAP, Feature plots  ##
#########################################
E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1 <- SCTransform(E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE, conserve.memory = TRUE)
E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1_pca <- RunPCA(E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1, verbose = TRUE)

pdf(file="./4-Subset1/ElbowPlot-PCA_E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1_pca.pdf", width=9, height=10)
ElbowPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1_pca, ndims = 50, reduction = "pca")
dev.off()
# 30 components

### Harmony integration by sequencing batch

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch <- RunHarmony(E875_SmoNull_Merge_RmDB_harmonySeqBatch__subset1_pca, 
				group.by.vars = "Library.ident", 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

# UMAP calculation and clustering
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch <- RunUMAP(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch, dims = 1:30, reduction = "harmony", verbose = TRUE)
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch <- FindNeighbors(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch,reduction = "harmony", dims = 1:30, verbose = TRUE)
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0 <- FindClusters(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch, verbose = TRUE, resolution = 1.0)

# UMAP DimPlots
pdf(file="./4-Subset1/UMAP_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T, raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T,  group.by = "Condition", raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, pt.size = 0.7, reduction = "umap", label = T,  group.by = "Sample", raster = FALSE))
dev.off()

# UMAP split by replicate
reps <- unique(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0$Replicate)
pdf(file = "./4-Subset1/UMAP-PerRep_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.pdf", width = 9, height = 10)
for (r in reps){
  print(DimPlot(subset(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, subset = Replicate == r), pt.size = 0.7, reduction = "umap", label = TRUE, group.by = "Condition", raster = FALSE) + 
  ggtitle(paste("Replicate:", r))
  )
  print(DimPlot(subset(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, subset = Replicate == r), pt.size = 0.7, reduction = "umap", label = TRUE, group.by = "Condition", raster = FALSE, cols = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C")) + 
  ggtitle(paste("Replicate:", r)) + 
  NoLegend()
  )
}
dev.off()

saveRDS(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.RDS")
#E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0 <- readRDS(file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.RDS")

########################################################
## Identification of the proper clustering resolution ##
########################################################
dataset <- readRDS(file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.RDS")

pdf(file="./4-Subset1/1-ResolutionLoop/E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_ResTest.pdf", width=9, height=10)
for(res0 in c(0.2,0.5,0.8,1,1.2,1.5,1.8,2,2.2,2.5,2.8,3)){
dataset_res0 <- FindClusters(dataset, verbose = TRUE, resolution = res0)
print(DimPlot(dataset_res0, raster=FALSE, reduction = "umap", label = T)+ggtitle(res0))
}
dev.off()

for(res in c(0.2,0.5,0.8,1,1.2,1.5,1.8,2,2.2,2.5,2.8,3)){
dataset_res <- FindClusters(dataset, verbose = TRUE, resolution = res)
# Find markers
dataset_res <- PrepSCTFindMarkers(dataset_res)
dataset_res.markers <- FindAllMarkers(dataset_res, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6,recorrect_umi = FALSE)
write.csv(dataset_res.markers, file = paste0("./4-Subset1/1-ResolutionLoop/E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res",res,"_all-markers.csv"))

# Marker heatmap
dataset_res.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = paste0("./4-Subset1/1-ResolutionLoop/TopMarkers-Heatmap_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res",res,".pdf"), width = 11, height = 11)
print(DoHeatmap(dataset_res, features = top15$gene) + NoLegend())
dev.off()

pdf(file = paste0("./4-Subset1/1-ResolutionLoop/TopMarkers-Dotplot_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res",res,".pdf"), width = 15, height = 40)
print(DotPlot(dataset_res, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

print(paste0("Done with resolution ",res))
}


# Interesting markers
gene_list <- c("Isl1","Osr1","Tbx1","Fgf10","Arg1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Nkx2-5","Tbx5","Gata4","Gli1","Gli2","Gli3","Foxf1","Hhip","Ptch1","Mab21l2","Hand1","Smarcd3","Shh","Ihh","Dhh","Alx1","Pecam1","Wt1","Tbx18")

pdf(file = "./4-Subset1/MarkerGenes-ORDER_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, raster=FALSE, features = i, order = TRUE) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()
pdf(file = "./4-Subset1/MarkerGenesByCond-ORDER_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res1.0, raster=FALSE, features = i, split.by = "Condition", keep.scale = "all", order = TRUE) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

## Continue with Resolution 2 ## 
	### Merge cluster ###
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data %>%
  mutate(Annotation = case_when(
     seurat_clusters %in% c(0,9,13,25,27,33,36,37) ~ "Mesenchyme",
	 seurat_clusters %in% c(1,2,3,4,11,16,17,19,21,22,28,32,35,39,40) ~ "Limb/Paraxial mesoderm",
	 seurat_clusters %in% 38 ~ "Neural crest",
	 seurat_clusters %in% 12 ~ "JCF-like",
	 seurat_clusters %in% c(7,4) ~ "pSHF",
	 seurat_clusters %in% 18 ~ "aSHF",
	 seurat_clusters %in% 20 ~ "Head muscle progenitor",
	 seurat_clusters %in% c(6,8,15,24,34) ~ "CM",
	 seurat_clusters %in% 10 ~ "Differentiating CM",
	 seurat_clusters %in% 5 ~ "Epicardium",
	 seurat_clusters %in% c(23,31) ~ "Intermediate mesoderm",
	 seurat_clusters %in% c(7,14,26,29,30) ~ "pSHF",
	TRUE ~  seurat_clusters
  ))

## Color codes ##
AnnotationNullColors = c('#a6b2ff','#d62728','#ff7f0e','#2ca02c','#1f77b4','#d7aa9e','#9467bd','#17becf','#e377c2','#8c564b','#8eba96')
names(AnnotationNullColors) = c("Limb/Paraxial mesoderm","CM","Differentiating CM","aSHF","pSHF","Epicardium","Neural crest","Intermediate mesoderm","JCF-like","Mesenchyme","Head muscle progenitor")

### Heatmap annotated UMAP
AnnoGenes <- c(
"Prrx1","Lef1","Twist1",#Limb/Paraxial mesoderm
"Tnnt2","Myh7","Ryr2","Ttn",#CMs
"Armh4","Mef2c",#Differentiating CMs
"Isl1","Tbx1","Fgf10",#aSHF
"Arg1","Osr1","Aldh1a2",#pSHF
"Foxf1","Hand1","Mab21l2",#JCF-like
"Col1a1","Col1a2","Dkk2","Pitx2",#Mesenchyme
"Msc","Lhx2","Kdr",#Head muscle progenitors
"Sox10","Tfap2a","Dlx2",#Neural crest
"Tbx18","Wt1",#Epicardium
"Pax2","Pax8","Lhx1"#Intermediate mesoderm
)

E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data$Annotation <- factor(
  E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data$Annotation,
  levels = c("Limb/Paraxial mesoderm","CM","Differentiating CM","aSHF","pSHF","JCF-like","Mesenchyme","Head muscle progenitor","Neural crest","Epicardium","Intermediate mesoderm")
)

E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- PrepSCTFindMarkers(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2)
expression_matrix <- GetAssayData(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, assay = "SCT", slot = "data")
annotation_labels <- setNames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data$Annotation, rownames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data))

averaged_expression <- expression_matrix[AnnoGenes, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression") %>%
  mutate(Annotation = annotation_labels[Cell]) %>% 
  filter(!is.na(Annotation)) %>%  # Remove rows with NA annotations
  group_by(Gene, Annotation) %>%
  summarize(MeanExpression = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Annotation, values_from = MeanExpression) %>%
  column_to_rownames(var = "Gene")

scaled_expression <- t(scale(t(as.matrix(averaged_expression)), center = TRUE, scale = TRUE))
scaled_expression <- scaled_expression[AnnoGenes, ]

bottom_annotation <- HeatmapAnnotation(
  Annotation = colnames(scaled_expression),
  col = list(Annotation = AnnotationNullColors),
  annotation_label = "Cell type",
  show_legend = FALSE 
)

pdf(file = "./4-Subset1/Heatmap_AveragedExpression_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.pdf", width = 10, height = 6.9)
ht <- Heatmap(
  matrix = as.matrix(scaled_expression),
  name = NULL,
  width = unit(4.6, "in"),
  col = colorRamp2(c(-1.5,0,2.5), hcl_palette = "viridis"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = NULL,
  row_title = NULL,
  column_names_rot = 45,
  bottom_annotation = bottom_annotation,
  heatmap_legend_param = list(
  title = "Relative expression \n (z-score)",
  at = c(-1.5, 0, 2.5),
  title_position = "leftcenter-rot",
  labels_gp = gpar(fontsize = 8),
  legend_direction = "vertical",
  title_gp = gpar(fontsize = 8, fontface = "bold")
  )
)
draw(ht, heatmap_legend_side = "left")
dev.off()

pdf(file="./4-Subset1/UMAP-Color_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "Annotation", label = T, cols = AnnotationNullColors))
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "Annotation", label = T, cols = AnnotationNullColors) + NoLegend())
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, pt.size = 0.7, reduction = "umap",  group.by = "Condition", cols = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C"), raster = FALSE))
print(DimPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, pt.size = 0.7, reduction = "umap",  group.by = "Condition", cols = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C"), raster = FALSE) + NoLegend())
dev.off()

# Add Log count per million normalization (RNA slot) for gene expression visualization and analysis
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- NormalizeData(
  E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1000000,
  verbose = TRUE
)

TmpObj_GeneExpr <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2
DefaultAssay(TmpObj_GeneExpr) <- "RNA"

GENES <- c("Mab21l2","Hand1","Foxf1",
		   "Tbx1","Fgf10","Mef2c",
		   "Isl1","Osr1","Aldh1a2",
		   "Ttn","Tnnt2","Ryr2",
		   "Ptch1","Hhip","Gli1")

pdf(file = "./4-Subset1/SupFigGENES-ORDER-logCPM_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.pdf", width = 10, height = 10)
for(i in GENES){			   
print(FeaturePlot(TmpObj_GeneExpr, raster=FALSE, features = i, order = TRUE) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

# Save
saveRDS(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.RDS")
#E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- readRDS(file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.RDS")

################################################################
## Cell composition analysis with Propeller (Speckle package) ## !! R/4.3.1 !!
################################################################
module load R/4.3.1

# R/4.3.1
library(Seurat)
library(speckle)
library(limma)
library(ggplot2)

E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0 <- readRDS(file = "./4-Subset1/3-CellProportion_analysis/E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0-WithAnnot.RDS")

propResults <- propeller(clusters=E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0$Subset1Annot, sample=E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0$Sample, group=E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0$Condition)
write.csv(propResults, file = "./4-Subset1/3-CellProportion_analysis/Propeller_FullDataset_AnnotationSubset1.csv")
#propResults <- read.csv("./4-Subset1/3-CellProportion_analysis/Propeller_FullDataset_AnnotationSubset1.csv")

## R/4.0.3 ##
propResults <- read.csv("./4-Subset1/3-CellProportion_analysis/Propeller_FullDataset_AnnotationSubset1.csv")

tab <- table(E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0@meta.data$Subset1Annot,E875_SmoNull_Merge_RmDB_harmonySeqBatch_res1.0@meta.data$Sample)
prop_tab <- prop.table(tab, margin = 1)
round(prop_tab, 3)

tab_df <- as.data.frame(tab) %>%
  rename(cluster = Var1, sample_id = Var2, count = Freq)
sample_totals <- colSums(tab)
tab_df <- tab_df %>%
  mutate(proportion = count / sample_totals[sample_id])

tab_df <- tab_df %>%
  mutate(group_id = case_when(
    sample_id %in% c("Sample1", "Sample3", "Sample5") ~ "Control",
    sample_id %in% c("Sample2", "Sample4", "Sample6") ~ "SmoNull"
  ))

celltypes <- c("pSHF", "aSHF", "JCF-like")
prep_subset <- tab_df %>% filter(cluster %in% celltypes)

prop_means <- propResults %>%
  filter(X %in% celltypes) %>%
  select(cluster = X, PropMean.Control, PropMean.SmoNull) %>%
  pivot_longer(cols = starts_with("PropMean"), names_to = "Condition", values_to = "PropMean") %>%
  mutate(Condition = ifelse(grepl("Control", Condition), "Control", "SmoNull"))

sd_df <- prep_subset %>%
  group_by(cluster, group_id) %>%
  summarise(SD = sd(proportion), .groups = "drop") %>%
  rename(Condition = group_id)
plot_df <- left_join(prop_means, sd_df, by = c("cluster", "Condition"))

pdf(file="./4-Subset1/3-CellProportion_analysis/Boxplot_JCF-SHF_proportions.pdf", width=10, height=8)
ggplot(prep_subset, aes(x = group_id, y = proportion)) +
  geom_boxplot(aes(color = group_id), fill = NA, width = 0.5, size = 1) +
  geom_jitter(color = "black", width = 0.1, size = 2) +
  scale_color_manual(values = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C")) +
  facet_wrap2(~ cluster, scales = "free_y", strip = strip_themed(
    background_x = elem_list_rect(fill = c("aSHF" = "#2ca02c", "pSHF" = "#1f77b4", "JCF-like" = "#e377c2")[levels(factor(prep_subset$cluster))], color = "gray"),
    text_x = elem_list_text(color = "white", face = "bold")
  )) +
  theme_light(base_size = 18) +
  theme(axis.ticks = element_blank(), legend.position = "none", panel.grid = element_blank(), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12), rect = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = "gray", linewidth = 1)) +
  ylab("Proportion") + xlab("")
dev.off()

###################################################
#########  Trajectory analysis slingshot  #########
###################################################
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- readRDS(file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.RDS")

### Subset for slingshot ###
E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM <- subset(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2,  subset =  Annotation %in% c("aSHF","Differentiating CM","CM"), invert = FALSE)

saveRDS(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.RDS")
#E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM <- readRDS(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.RDS")

pdf(file="./4-Subset1/4-Slingshot/1-Subset_aSHF/UMAP-Color_E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "Annotation", label = T, cols = AnnotationNullColors))
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "Annotation", label = T, cols = AnnotationNullColors) + NoLegend())
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend())
dev.off()

# Figure
E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data <- E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data %>%
  mutate(AnnotationCM = case_when(
     seurat_clusters %in% 15 ~ "Outflow tract CM",
     seurat_clusters %in% 6 ~ "Ventricular CM", 
	 seurat_clusters %in% c(8,34) ~ "Atrial CM",
	 seurat_clusters %in% 24 ~ "Venous pole",
	 seurat_clusters %in% 18 ~ "aSHF",
	 seurat_clusters %in% 10 ~ "Differentiating CM",
	TRUE ~  seurat_clusters
  ))

## Color codes ##
AnnotationCMColors = c('#ff7f0e','#2ca02c','#EB531B','#D62728','#A12F3C','#6B3750')
names(AnnotationCMColors) = c("Differentiating CM","aSHF","Outflow tract CM","Ventricular CM","Atrial CM","Venous pole")

pdf(file="./4-Subset1/4-Slingshot/1-Subset_aSHF/UMAP-Figure_E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.pdf", width=9, height=10)
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "AnnotationCM", label = T, cols = AnnotationCMColors))
print(DimPlot(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, pt.size = 0.7, raster=FALSE, reduction = "umap", group.by = "AnnotationCM", label = T, cols = AnnotationCMColors) + NoLegend())
dev.off()

### Heatmap annotated UMAP
AnnoGenes <- c(
"Tbx1","Fgf10","Isl1",#aSHF
"Armh4",#Differentiating CMs
"Tdgf1", #OFT
"Irx4","Nppa","Bmp10", # Ventricular CM
"Tbx2","Nr2f2", # Atrial CM
"Wnt2","Hcn4","Shox2" # Venous pole
)

E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$AnnotationCM <- factor(
  E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$AnnotationCM,
  levels = c("aSHF","Differentiating CM","Outflow tract CM","Ventricular CM","Atrial CM","Venous pole")
)

E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM <- PrepSCTFindMarkers(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM)
expression_matrix <- GetAssayData(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM, assay = "SCT", slot = "data")
annotation_labels <- setNames(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$AnnotationCM, rownames(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data))

averaged_expression <- expression_matrix[AnnoGenes, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Cell", values_to = "Expression") %>%
  mutate(AnnotationCM = annotation_labels[Cell]) %>% 
  filter(!is.na(AnnotationCM)) %>%  # Remove rows with NA annotations
  group_by(Gene, AnnotationCM) %>%
  summarize(MeanExpression = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = AnnotationCM, values_from = MeanExpression) %>%
  column_to_rownames(var = "Gene")

scaled_expression <- t(scale(t(as.matrix(averaged_expression)), center = TRUE, scale = TRUE))
scaled_expression <- scaled_expression[AnnoGenes, ]

bottom_annotation <- HeatmapAnnotation(
  AnnotationCM = colnames(scaled_expression),
  col = list(AnnotationCM = AnnotationCMColors),
  annotation_label = "Cell type",
  show_legend = FALSE 
)

pdf(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/Heatmap_AveragedExpression_E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.pdf", width = 10, height = 6.9)
ht <- Heatmap(
  matrix = as.matrix(scaled_expression),
  name = NULL,
  width = unit(4.6, "in"),
  col = colorRamp2(c(-1.5,0,2.5), hcl_palette = "viridis"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = NULL,
  row_title = NULL,
  column_names_rot = 45,
  bottom_annotation = bottom_annotation,
  heatmap_legend_param = list(
  title = "Relative expression \n (z-score)",
  at = c(-1.5, 0, 2.5),
  title_position = "leftcenter-rot",
  labels_gp = gpar(fontsize = 8),
  legend_direction = "vertical",
  title_gp = gpar(fontsize = 8, fontface = "bold")
  )
)
draw(ht, heatmap_legend_side = "left")
dev.off()

##START c18 END c6 HARMONY
sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30 <- as.SingleCellExperiment(E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM)
sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30 <- slingshot(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30, clusterLabels = E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM$seurat_clusters, reducedDim = reducedDim(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30, "HARMONY")[,1:30], approx_points = 200, start.clus="18", end.clus="6", extend = 'n', stretch = 0)

saveRDS(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30, file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30-extdN-Stretch0.RDS")
#sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30 <- readRDS(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30-extdN-Stretch0.RDS")

#Check the identified lineages and Curves
SlingshotDataSet(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)

colnames(colData(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30))

	# Plotting pseudotimes
colors <- colorRampPalette(c("#0D0887FF", "#CC4678FF", "#F0F921FF"))(100)
pdf(file="./4-Subset1/4-Slingshot/1-Subset_aSHF/sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30_One-functionMethod_Pseudotimes-extdN-Stretch0nolinesFIGURE.pdf", width=12, height=12)
for(u in seq(24,26)){
plotcol <- colors[cut(as.matrix(colData(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)[u]), breaks=100)]
plotcol[is.na(plotcol)] <- "#D9D9D9"
plot(reducedDims(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)$UMAP, col = plotcol, pch=16, asp = 1, main=paste0(colnames(colData(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)[u])))
#lines(SlingshotDataSet(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30), lwd=2, col='black')
}
dev.off()
# lineage 3 is the good one

#==# Selected cells only for pseudotime #==#
GoodPseudoCells <- rownames(na.omit(colData(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)["slingPseudotime_3"]))

E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2[["CellName"]] <- colnames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2)
aSHFTraj_Slingshot_SelectedCells <- subset(x = E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, subset = CellName %in% GoodPseudoCells, invert=FALSE)

meta <- na.omit(colData(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30)[,c("seurat_clusters", "Condition", "Annotation", "slingPseudotime_3")])
colnames(meta) <- c("seurat_clusters", "Condition", "Annotation", "dptval")
aSHFTraj_Slingshot_SelectedCells[["SlingshotPseudotime"]] <- meta$dptval[match(rownames(aSHFTraj_Slingshot_SelectedCells@meta.data), rownames(meta))]

### Save seurat object with only selected cells ###
saveRDS(aSHFTraj_Slingshot_SelectedCells, file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/aSHFTraj_Slingshot_SelectedCells.RDS")
#aSHFTraj_Slingshot_SelectedCells <- readRDS(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/aSHFTraj_Slingshot_SelectedCells.RDS")
######

pdf(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/TrajectoryViz_aSHF.pdf", width = 10, height = 10)
  good2 <- data.frame(cells = colnames(aSHFTraj_Slingshot_SelectedCells), Pseudotime = aSHFTraj_Slingshot_SelectedCells@meta.data$SlingshotPseudotime)
  all2 <- data.frame(cells = rownames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data), Pseudotime = NA)
  rest2 <- all2[!(all2$cells %in% good2$cells),]
  full2 <- rbind(good2,rest2)
temp <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2
temp@meta.data$Pseudo <- full2$Pseudotime[match(rownames(temp@meta.data), full2$cells)]
print(FeaturePlot(temp, features = "Pseudo", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
dev.off()

## Differential Progression ##
library(patchwork)

psts <- slingPseudotime(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$Condition,
		 annotations = E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$Annotation) %>%
  pivot_longer(starts_with("curve"), values_to = "pseudotime", names_to = "lineages")

pdf(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/DiffProgression_visualisation.pdf", width = 10, height = 10)
print(ggplot(psts, aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")
  )
print(ggplot(psts[psts$lineages == "curve3",], aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")
  )
p_density <- ggplot(psts[psts$lineages == "curve3",], aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C")) +
  facet_wrap(~lineages) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank()
  )
p_anno <- ggplot(psts[psts$lineages == "curve3",], aes(x = pseudotime, y = 1, fill = annotations)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "CM" = "#d62728",
    "Differentiating CM" = "#ff7f0e",
    "aSHF" = "#2ca02c"
  )) +
  theme_void() +
  theme(legend.position = "bottom")
combined_plot <- p_anno / p_density + plot_layout(heights = c(0.2, 1))
print(combined_plot)
dev.off()


	# Condiment test # R/4.3.1 #
module load ggc
module load python/anaconda-2022.05
module load R/4.3.1

  library(Seurat)
  library(ggplot2)
  library(slingshot)
  library(condiments)

E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM <- readRDS(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM.RDS")
sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30 <- readRDS(file = "./4-Subset1/4-Slingshot/1-Subset_aSHF/sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30-extdN-Stretch0.RDS")
  
prog_res <- progressionTest(sce_SmoNull_subset.aSHF.CM_startc18endc6_HARMONY30, conditions = E875_SmoNull_Merge_RmDB_harmonySeqBatch_subset.aSHF.CM@meta.data$Condition, global = FALSE, lineages = TRUE)
knitr::kable(prog_res)

#|lineage | statistic| p.value|
#|:-------|---------:|-------:|
#|1       | 0.3106866|       0|
#|2       | 0.2277300|       0|
#|3       | 0.1303265|       0|


######################################
## UCell function Hh score # R/4.3.1 #
######################################
# R/4.3.1
library(Seurat)
library(UCell)
library(ggplot2)
library(circlize)
library(ggh4x)

E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- readRDS(file = "./E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.RDS")
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2 <- PrepSCTFindMarkers(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2)

markers <- list(Limited_Hh_list = c("Ptch1","Hhip","Gli1"))

SeuratObj <- AddModuleScore_UCell(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, markers, assay = "SCT")
SeuratObj <- SmoothKNN(SeuratObj, signature.names = paste0(names(markers),"_UCell"), reduction="harmony", k=10)

UCell_score <- data.frame(UCell = SeuratObj@meta.data$Limited_Hh_list_UCell_kNN, row.names = rownames(SeuratObj@meta.data))
write.csv(UCell_score, file = "./4-Subset1/5-UCell/Limited_Hh_list_UCell_kNN_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2_AfterPrepSCTFindMarkers.csv")

pdf(file="./4-Subset1/5-UCell/UMAP_UCell-SCT_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2_AfterPrepSCTFindMarkers.pdf", width=9, height=10)
color_scale <- colorRamp2(c(0, 0.02 ,0.05), hcl_palette = "cividis")
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2, "cm"))
print(FeaturePlot(SeuratObj,features = paste0(names(markers),"_UCell"), order = TRUE, coord.fixed = F, pt.size = 0.7) +
  theme(axis.line = element_line(arrow = arrow(length=unit(0.1, "inches"))), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(hjust = 0)) +
  labs(x="UMAP1", y="UMAP2", title= paste0(names(markers),"_UCell")) +
  guides(x = axis, y= axis) +
  scale_color_gradientn(colors = color_scale(seq(0, 0.1, length.out = 10)), breaks = scales::breaks_pretty(n = 3))
)
dev.off()

################## R/4.0.3 ##########
UCell_score <- read.csv("./4-Subset1/5-UCell/Limited_Hh_list_UCell_kNN_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2_AfterPrepSCTFindMarkers.csv")

E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2[["Limited_Hh_list_UCell"]] <- UCell_score$UCell[match(rownames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data), UCell_score$X)]

pdf(file = "./4-Subset1/5-UCell/HhScore-UMAP_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.pdf", width = 10, height = 10)	   
print(FeaturePlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, features = "Limited_Hh_list_UCell", order = TRUE, pt.size = 0.6) + 
	  scale_color_gradientn(colors = c("#512939", "#F38D68", "#F38D68"),
      values = scales::rescale(c(0, 0.4, 1)),
      limits = c(0, 1),
      na.value = "#E4E4E4"
    )
)
dev.off()	

### Hh score Vlnplot per cell type
UCell_score <- read.csv("./4-Subset1/5-UCell/Limited_Hh_list_UCell_kNN_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2_AfterPrepSCTFindMarkers.csv")
E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2[["Limited_Hh_list_UCell"]] <- UCell_score$UCell[match(rownames(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data), UCell_score$X)]

pdf(file = "./4-Subset1/5-UCell/ViolinPlot_HhScore-perCellType_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.pdf", width = 11, height = 11)	   
Idents(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2) <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data$Condition
print(VlnPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, "Limited_Hh_list_UCell", pt.size = 0, split.by = "Condition")+ 
	    NoLegend() +
		scale_fill_manual(values = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C")) +
		stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.2, color = "white", size = 0.2) +
		labs(title = NULL,
			     x = NULL,              
                 y = "Hedgehog Score (AU)") + 
		theme(axis.title.x = element_blank(),
		axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)))
Idents(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2) <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data$Annotation
print(VlnPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, "Limited_Hh_list_UCell", pt.size = 0)+ 
	    NoLegend() +
		geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
		scale_fill_manual(values = AnnotationNullColors) +
		labs(title = NULL,
			     x = NULL,              
                 y = "Hedgehog Score (AU)") + 
		theme(axis.title.x = element_blank(),
		axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)))
print(VlnPlot(E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2, "Limited_Hh_list_UCell", pt.size = 0, split.by = "Condition")+ 
	    NoLegend() +
		scale_fill_manual(values = c("Control" = "#D3D3D3", "SmoNull" = "#DE125C")) +
		stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.2, color = "white", size = 0.2) +
		labs(title = NULL,
			     x = NULL,              
                 y = "Hedgehog Score (AU)") + 
		theme(axis.title.x = element_blank(),
		axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)))
dev.off()

library(rstatix)

# Test for global score Ctl vs. SmoNull
WT <- wilcox_test(data = E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data, formula = Limited_Hh_list_UCell ~ Condition, p.adjust.method = "fdr")
write.csv(WT, file = "./4-Subset1/5-UCell/Wilcox-GlobalPerCond_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.csv")
#.y.                   group1  group2     n1    n2 statistic     p
#* <chr>                 <chr>   <chr>   <int> <int>     <dbl> <dbl>
#1 Limited_Hh_list_UCell Control SmoNull 13620 11998 131705389     0

# Test Ctl vs. SmoNull per cell type
PWT <- E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2@meta.data %>%
  group_by(Annotation) %>%
  filter(n_distinct(Condition) == 2) %>%  # Only keep groups with both conditions
  wilcox_test(Limited_Hh_list_UCell ~ Condition) %>% 
  adjust_pvalue(method='fdr')
  
write.csv(PWT, file = "./4-Subset1/5-UCell/ParwiseWilcox-PerCluster_E875_SmoNull_Merge_RmDB_subset1_harmonySeqBatch_res2.csv")

