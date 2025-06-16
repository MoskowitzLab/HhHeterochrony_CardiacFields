########################################################################################################################################################################
################################################                                                                          ##############################################
################################################ Krup-Pijuan-Tyser integrated Atlas - Tyser pre-processing (v03-31-2025)  ##############################################
################################################                                                                          ##############################################
########################################################################################################################################################################
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

#======================================================================================================================================================================#
######################################################========       Data loading and preparation        ========#######################################################
#======================================================================================================================================================================#
# https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/

sceNormCluster <- readRDS(file = "./sce_goodQual.NORM.clusters.Rds")

# Conversion to Seurat object
seurat_object <- CreateSeuratObject(
    counts = counts(sceNormCluster),
    meta.data = as.data.frame(colData(sceNormCluster)),
    project = "MyProject"
)

# Add logcounts as an assay to the Seurat object
seurat_object[["RNA"]] <- CreateAssayObject(counts = counts(sceNormCluster))
seurat_object[["RNA"]]@data <- logcounts(sceNormCluster)

# Convert UMAP embeddings to matrix and ensure proper column naming
if ("UMAP" %in% reducedDimNames(sceNormCluster)) {
    umap_data <- as.matrix(reducedDim(sceNormCluster, "UMAP")[,-3])
    colnames(umap_data) <- paste0("UMAP_", 1:ncol(umap_data))
    seurat_object[["umap"]] <- CreateDimReducObject(
        embeddings = umap_data,
        key = "UMAP_",
        assay = "RNA"
    )
}

# Set the active identities to the clusters
Idents(seurat_object) <- seurat_object@meta.data$cluster

# Version data save
saveRDS(seurat_object, file = "./seurat_object_goodQual.NORM.clusters.RDS")
#seurat_object <- readRDS(file = "./seurat_object_goodQual.NORM.clusters.RDS")

#======================================================================================================================================================================#
################################################========                Subset relevant cells                 ========##################################################
#======================================================================================================================================================================#
seurat_object_goodQual <- readRDS(file = "./seurat_object_goodQual.NORM.clusters.RDS")
sceTyserRaw <- readRDS(file = "./heartData_unbiased.RAW.Rds") # raw data matrix

# Select only good quality cells from Tyser paper analysis
goodCells <- colnames(seurat_object_goodQual)
sceTyserRaw_filt1 <- sceTyserRaw[,colnames(sceTyserRaw) %in% goodCells]

# Remove unwanted gene info + ERCC spike-in rows
remove_rows <- c('__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique')
sceTyserRaw_filt1 <- sceTyserRaw_filt1[!(rownames(sceTyserRaw_filt1) %in% remove_rows | grepl('ERCC-', rownames(sceTyserRaw_filt1))), ]

# Convert ENSMUS gene codes to Symbols
library("AnnotationDbi")
library("org.Mm.eg.db")

Genes=as.data.frame(rownames(sceTyserRaw_filt1))
colnames(Genes)<-"ensid"
Genes$symb = mapIds(org.Mm.eg.db,
                    keys=Genes$ensid, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

	# Dealing with duplicates
dup=table(Genes$symb)
dup2 <- dup[dup > 1]
dim(dup2) # 39 duplicates

DupInfo <- Genes[Genes$symb %in% as.data.frame(dup2)$Var1,]

# Merge the expression matrix with the gene table based on EnsemblID
sceTyserRaw_filt2 <- merge(as.data.frame(sceTyserRaw_filt1), Genes, by.x = 'row.names', by.y = 'ensid')
rownames(sceTyserRaw_filt2) <- sceTyserRaw_filt2$Row.names
sceTyserRaw_filt2$Row.names <- NULL

# Sum expression values by gene symbol
summed_expression <- aggregate(. ~ symb, data = sceTyserRaw_filt2, FUN = sum)

# Convert back to matrix format
sceTyserRaw_filt2_mat <- as.matrix(summed_expression[, -1])
rownames(sceTyserRaw_filt2_mat) <- summed_expression$symb

# Create subset and curated Seurat object
SeuratTyser <- CreateSeuratObject(counts = sceTyserRaw_filt2_mat, meta.data = as.data.frame(seurat_object_goodQual@meta.data), min.cells = 3, min.features = 200)

# Version data save
saveRDS(SeuratTyser, file = "./Integration/SeuratTyser.RDS")
#SeuratTyser <- readRDS(file = "./Integration/SeuratTyser.RDS")

#########################
## Cell cycle scoring  ##
#########################
cc.genes <- read.table("./0-R_codes/2-Ressources/Cell-Cycle_genes_scRNA-Seq_MusMusculus.txt", sep="\t", header = T,fill = TRUE)
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[cc.genes$CellCyclePhase == "S",1]
g2m.genes <- cc.genes[cc.genes$CellCyclePhase == "G2/M",1]

# Assign cell cycle scores
SeuratTyser <- SCTransform(SeuratTyser, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = TRUE)
SeuratTyser <- CellCycleScoring(SeuratTyser, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = SeuratTyser@meta.data)

saveRDS(SeuratTyser, file = "./Integration/SeuratTyser_CCscoring.RDS")
#SeuratTyser <- readRDS(file = "./Integration/SeuratTyser_CCscoring.RDS")

#############################################
## Final SCT - Cleaned individual datasets ##
#############################################
SeuratTyser <- SCTransform(SeuratTyser, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)

saveRDS(SeuratTyser, file = "./Integration/SeuratTyser_SCT.RDS")
#SeuratTyser <- readRDS(file = "./Integration/SeuratTyser_SCT.RDS")

#======================================================================================================================================================================#
########################################========            Subset 1 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
SeuratTyser <- readRDS(file = "./Integration/SeuratTyser_CCscoring.RDS")
CellsToRemove_Tyser <- read.csv("./Krup-Tyser-Pijuan_Joint/1-Merge1/PrepSub1_CellsToRemove_Tyser.csv")

SeuratTyser_subForMerge_Krup_Tyser_Pijuan <- subset(SeuratTyser, cells = CellsToRemove_Tyser$x, invert=TRUE)
SeuratTyser_subForMerge_Krup_Tyser_Pijuan <- SCTransform(SeuratTyser_subForMerge_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(SeuratTyser_subForMerge_Krup_Tyser_Pijuan, file = "./Integration/SeuratTyser_subForMerge_Krup_Tyser_Pijuan.RDS")
#SeuratTyser_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./Integration/SeuratTyser_subForMerge_Krup_Tyser_Pijuan.RDS")

#======================================================================================================================================================================#
########################################========            Subset 2 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
	# In addition to previous 'subset1'
SeuratTyser_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./Integration/SeuratTyser_subForMerge_Krup_Tyser_Pijuan.RDS")
CellsToRemove2_Tyser <- read.csv("./Krup-Tyser-Pijuan_Joint/2-Merge2/PrepSub2_CellsToRemove_Tyser.csv")

SeuratTyser_subForMerge3_Krup_Tyser_Pijuan <- subset(SeuratTyser_subForMerge_Krup_Tyser_Pijuan, cells = CellsToRemove2_Tyser$x, invert=TRUE)
SeuratTyser_subForMerge3_Krup_Tyser_Pijuan <- SCTransform(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(SeuratTyser_subForMerge3_Krup_Tyser_Pijuan, file = "./Integration/SeuratTyser_subForMerge3_Krup_Tyser_Pijuan.RDS")
#SeuratTyser_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./Integration/SeuratTyser_subForMerge3_Krup_Tyser_Pijuan.RDS")

#43076 features across 2077 samples within 2 assays
#Active assay: SCT (20514 features, 3000 variable features)
