########################################################################################################################################################################
################################################                                                                          ##############################################
################################################ Krup-Pijuan-Tyser integrated Atlas - Pijuan pre-processing (v03-31-2025) ##############################################
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
  
  library(MouseGastrulationData)

#======================================================================================================================================================================#
######################################################========       Data loading and preparation        ========#######################################################
#======================================================================================================================================================================#

# Loading of the complete WT dataset
sce <- EmbryoAtlasData(type = "processed")
saveRDS(sce, file = "./sce_object.RDS")
#sce <- readRDS(file = "./sce_object.RDS")

# Add log-norm counts
library(scuttle)
sce2 <- logNormCounts(sce)

# Conversion to Seurat object
seurat_object <- CreateSeuratObject(
    counts = counts(sce2),
    meta.data = as.data.frame(colData(sce2)),
    project = "MyProject"
)

# Add logcounts as an assay to the Seurat object
seurat_object[["RNA"]] <- CreateAssayObject(counts = counts(sce2))
seurat_object[["RNA"]]@data <- as.matrix(logcounts(sce2))

# Add UMAP embeddings to matrix and ensure proper column naming
if ("umap" %in% reducedDimNames(sce2)) {
    umap_data <- as.matrix(reducedDim(sce2, "umap")[,-3])
    colnames(umap_data) <- paste0("UMAP_", 1:ncol(umap_data))
    seurat_object[["umap"]] <- CreateDimReducObject(
        embeddings = umap_data,
        key = "UMAP_",
        assay = "RNA"
    )
}

# Set the active identities to the clusters
Idents(seurat_object) <- seurat_object@meta.data$celltype

# Version data save
saveRDS(seurat_object, file = "./seurat_object_Gastrulation_Processed_singlets.celltype.RDS")
#seurat_object <- readRDS(file = "./seurat_object_Gastrulation_Processed_singlets.celltype.RDS")

#======================================================================================================================================================================#
################################################========                Subset relevant cells                 ========##################################################
#======================================================================================================================================================================#
Subseurat_object <- subset(seurat_object, idents = c("Epiblast","Primitive Streak","Anterior Primitive Streak","Nascent mesoderm","Mixed mesoderm","ExE mesoderm","Mesenchyme","Pharyngeal mesoderm","Cardiomyocytes"), invert=FALSE)
Subseurat_object@active.ident <- factor(x = Subseurat_object@active.ident, levels = c("Epiblast","Primitive Streak","Anterior Primitive Streak","Nascent mesoderm","Mixed mesoderm","ExE mesoderm","Mesenchyme","Pharyngeal mesoderm","Cardiomyocytes"))

Metadata_GoodQualMesoCells <- Subseurat_object@meta.data[Subseurat_object@meta.data$celltype %in% c("Primitive Streak","Anterior Primitive Streak","Nascent mesoderm","Mixed mesoderm","ExE mesoderm","Mesenchyme","Pharyngeal mesoderm","Cardiomyocytes"),]

# Select only good quality cells from the original paper analysis
Counts_GoodQualMesoCells <- assay(sce,"counts")[,colnames(assay(sce,"counts")) %in% rownames(Metadata_GoodQualMesoCells)]

# Convert ENSMUS genes to Symbols
library("AnnotationDbi")
library("org.Mm.eg.db")

Genes=as.data.frame(rownames(Counts_GoodQualMesoCells))
colnames(Genes)<-"ensid"
Genes$symb = mapIds(org.Mm.eg.db,
                    keys=Genes$ensid, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")			
	# Check for NA values in the symb column
sum(is.na(Genes$symb))

	# Dealing with duplicates
dup=table(Genes$symb)
dup2 <- dup[dup > 1]
dim(dup2) # 37 duplicates

DupInfo <- Genes[Genes$symb %in% as.data.frame(dup2)$Var1,]

# Merge the expression matrix with the gene table based on EnsemblID
Counts_GoodQualMesoCells_filt <- merge(as.data.frame(Counts_GoodQualMesoCells), Genes, by.x = 'row.names', by.y = 'ensid')
rownames(Counts_GoodQualMesoCells_filt) <- Counts_GoodQualMesoCells_filt$Row.names
Counts_GoodQualMesoCells_filt$Row.names <- NULL

# Sum expression values by gene symbol
summed_expression <- aggregate(. ~ symb, data = Counts_GoodQualMesoCells_filt, FUN = sum)

# Convert back to matrix format
GastrulationMeso_mat <- as.matrix(summed_expression[, -1])
rownames(GastrulationMeso_mat) <- summed_expression$symb

# Create subset and curated Seurat object
SeuratGastrulationMeso <- CreateSeuratObject(counts = GastrulationMeso_mat, meta.data = as.data.frame(Metadata_GoodQualMesoCells), min.cells = 3, min.features = 200)

# Remove cells from 'Mixed gastrulation'
SeuratGastrulationMeso2 <- subset(SeuratGastrulationMeso, subset = stage == "mixed_gastrulation", invert=TRUE)
# 24,921 cells

#########################
## Cell cycle scoring  ##
#########################
cc.genes <- read.table("./0-R_codes/2-Ressources/Cell-Cycle_genes_scRNA-Seq_MusMusculus.txt", sep="\t", header = T,fill = TRUE)
	# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[cc.genes$CellCyclePhase == "S",1]
g2m.genes <- cc.genes[cc.genes$CellCyclePhase == "G2/M",1]

# Assign cell cycle scores
SeuratGastrulationMeso2 <- SCTransform(SeuratGastrulationMeso2, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = TRUE)
SeuratGastrulationMeso2 <- CellCycleScoring(SeuratGastrulationMeso2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = SeuratGastrulationMeso2@meta.data)

## change orig.ident by individual sample names
SeuratGastrulationMeso2@meta.data$orig.ident <- SeuratGastrulationMeso2@meta.data$sample

#############################################
## Final SCT - Cleaned individual datasets ##
#############################################
SeuratGastrulationMeso2 <- SCTransform(SeuratGastrulationMeso2, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)
saveRDS(SeuratGastrulationMeso2, file = "./SeuratGastrulationMeso_SCT.RDS")
#SeuratGastrulationMeso_SCT <- readRDS(file = "./SeuratGastrulationMeso_SCT.RDS")


#======================================================================================================================================================================#
########################################========            Subset 1 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
SeuratGastrulationMeso2 <- readRDS(file = "./SeuratGastrulationMeso2.RDS")
CellsToRemove_Pijuan <- read.csv("./Krup-Tyser-Pijuan_Joint/1-Merge1/PrepSub1_CellsToRemove_Pijuan.csv")

SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan <- subset(SeuratGastrulationMeso2, cells = CellsToRemove_Pijuan$x, invert=TRUE)
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan <- SCTransform(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan, file = "./SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan.RDS")
#SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan.RDS")

#======================================================================================================================================================================#
########################################========            Subset 2 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
	# In addition to previous 'subset1'
SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan.RDS")
CellsToRemove2_Pijuan <- read.csv("./Krup-Tyser-Pijuan_Joint/1-Merge1/PrepSub1_CellsToRemove_Pijuan.csv")

SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan <- subset(SeuratGastrulationMeso2_subForMerge_Krup_Tyser_Pijuan, cells = CellsToRemove2_Pijuan$x, invert=TRUE)
SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan <- SCTransform(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan, file = "./SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan.RDS")
#SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./SeuratGastrulationMeso2_subForMerge3_Krup_Tyser_Pijuan.RDS")

#34037 features across 21565 samples within 2 assays
#Active assay: SCT (16636 features, 3000 variable features)


















