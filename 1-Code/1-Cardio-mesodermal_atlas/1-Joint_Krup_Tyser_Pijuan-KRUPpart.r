########################################################################################################################################################################
################################################                                                                        ################################################
################################################ Krup-Pijuan-Tyser integrated Atlas - Krup pre-processing (v03-14-2025) ################################################
################################################                                                                        ################################################
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
######################################################========      Data preparation ## v1.1.0-Prep      ========#######################################################
#======================================================================================================================================================================#
# sample ALK05_E7.75_rep1 removed for low quality

# Load the different datasets
ALK06_E60_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK06_E60_rep1/GSM6338196_ALK06_2_E60_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK06_E60_rep1/GSM6338196_ALK06_2_E60_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK06_E60_rep1/GSM6338196_ALK06_2_E60_con_rep1_barcodes.tsv.gz")
ALK06_E60_rep2.data <- ReadMtx(mtx = "./1-RawData/ALK06_E60_rep2/GSM6338197_ALK06_4_E60_con_rep2_matrix.mtx.gz", features = "./1-RawData/ALK06_E60_rep2/GSM6338197_ALK06_4_E60_con_rep2_features.tsv.gz", cells = "./1-RawData/ALK06_E60_rep2/GSM6338197_ALK06_4_E60_con_rep2_barcodes.tsv.gz")
ALK08_E60_rep3.data <- ReadMtx(mtx = "./1-RawData/ALK08_E60_rep3/GSM6338200_ALK08_20_E60_con_rep3_matrix.mtx.gz", features = "./1-RawData/ALK08_E60_rep3/GSM6338200_ALK08_20_E60_con_rep3_features.tsv.gz", cells = "./1-RawData/ALK08_E60_rep3/GSM6338200_ALK08_20_E60_con_rep3_barcodes.tsv.gz")
ALK08_lateE60_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK08_lateE60_rep1/GSM6338201_ALK08_14_lateE60_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK08_lateE60_rep1/GSM6338201_ALK08_14_lateE60_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK08_lateE60_rep1/GSM6338201_ALK08_14_lateE60_con_rep1_barcodes.tsv.gz")
ALK07_E65_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK07_E65_rep1/GSM6338198_ALK07_15_E65_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK07_E65_rep1/GSM6338198_ALK07_15_E65_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK07_E65_rep1/GSM6338198_ALK07_15_E65_con_rep1_barcodes.tsv.gz")
ALK07_lateE65_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK07_lateE65_rep1/GSM6338203_ALK07_3_lateE65_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK07_lateE65_rep1/GSM6338203_ALK07_3_lateE65_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK07_lateE65_rep1/GSM6338203_ALK07_3_lateE65_con_rep1_barcodes.tsv.gz")
ALK07_E70_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK07_E70_rep1/GSM6338204_ALK07_14_E70_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK07_E70_rep1/GSM6338204_ALK07_14_E70_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK07_E70_rep1/GSM6338204_ALK07_14_E70_con_rep1_barcodes.tsv.gz")
ALK07_earlyE75_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK07_earlyE75_rep1/GSM6338205_ALK07_7_earlyE75_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK07_earlyE75_rep1/GSM6338205_ALK07_7_earlyE75_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK07_earlyE75_rep1/GSM6338205_ALK07_7_earlyE75_con_rep1_barcodes.tsv.gz")
ALK07_lateE75_rep1.data <- ReadMtx(mtx = "./1-RawData/ALK07_lateE75_rep1/GSM6338207_ALK07_6_lateE75_con_rep1_matrix.mtx.gz", features = "./1-RawData/ALK07_lateE75_rep1/GSM6338207_ALK07_6_lateE75_con_rep1_features.tsv.gz", cells = "./1-RawData/ALK07_lateE75_rep1/GSM6338207_ALK07_6_lateE75_con_rep1_barcodes.tsv.gz")
ALK04_lateE75_rep2.data <- ReadMtx(mtx = "./1-RawData/ALK04_lateE75_rep2/GSM6338208_ALK04_3_lateE75_con_rep2_matrix.mtx.gz", features = "./1-RawData/ALK04_lateE75_rep2/GSM6338208_ALK04_3_lateE75_con_rep2_features.tsv.gz", cells = "./1-RawData/ALK04_lateE75_rep2/GSM6338208_ALK04_3_lateE75_con_rep2_barcodes.tsv.gz")
ALK07_E775_rep2.data <- ReadMtx(mtx = "./1-RawData/ALK07_E775_rep2/GSM6338210_ALK07_8_E775_con_rep2_matrix.mtx.gz", features = "./1-RawData/ALK07_E775_rep2/GSM6338210_ALK07_8_E775_con_rep2_features.tsv.gz", cells = "./1-RawData/ALK07_E775_rep2/GSM6338210_ALK07_8_E775_con_rep2_barcodes.tsv.gz")

# Create all Seurat objects (NOT AN ISSUE -> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-'))
ALK06_E60_rep1 <- CreateSeuratObject(counts = ALK06_E60_rep1.data, project = "ALK06_E60_rep1", min.cells = 3, min.features = 200)
ALK06_E60_rep2 <- CreateSeuratObject(counts = ALK06_E60_rep2.data, project = "ALK06_E60_rep2", min.cells = 3, min.features = 200)
ALK08_E60_rep3 <- CreateSeuratObject(counts = ALK08_E60_rep3.data, project = "ALK08_E60_rep3", min.cells = 3, min.features = 200)
ALK08_lateE60_rep1 <- CreateSeuratObject(counts = ALK08_lateE60_rep1.data, project = "ALK08_lateE60_rep1", min.cells = 3, min.features = 200)
ALK07_E65_rep1 <- CreateSeuratObject(counts = ALK07_E65_rep1.data, project = "ALK07_E65_rep1", min.cells = 3, min.features = 200)
ALK07_lateE65_rep1 <- CreateSeuratObject(counts = ALK07_lateE65_rep1.data, project = "ALK07_lateE65_rep1", min.cells = 3, min.features = 200)
ALK07_E70_rep1 <- CreateSeuratObject(counts = ALK07_E70_rep1.data, project = "ALK07_E70_rep1", min.cells = 3, min.features = 200)
ALK07_earlyE75_rep1 <- CreateSeuratObject(counts = ALK07_earlyE75_rep1.data, project = "ALK07_earlyE75_rep1", min.cells = 3, min.features = 200)
ALK07_lateE75_rep1 <- CreateSeuratObject(counts = ALK07_lateE75_rep1.data, project = "ALK07_lateE75_rep1", min.cells = 3, min.features = 200)
ALK04_lateE75_rep2 <- CreateSeuratObject(counts = ALK04_lateE75_rep2.data, project = "ALK04_lateE75_rep2", min.cells = 3, min.features = 200)
ALK07_E775_rep2 <- CreateSeuratObject(counts = ALK07_E775_rep2.data, project = "ALK07_E775_rep2", min.cells = 3, min.features = 200)

# Rename cells
ALK06_E60_rep1 <- RenameCells(object = ALK06_E60_rep1, add.cell.id = "ALK06_E60_rep1")
ALK06_E60_rep2 <- RenameCells(object = ALK06_E60_rep2, add.cell.id = "ALK06_E60_rep2")
ALK08_E60_rep3 <- RenameCells(object = ALK08_E60_rep3, add.cell.id = "ALK08_E60_rep3")
ALK08_lateE60_rep1 <- RenameCells(object = ALK08_lateE60_rep1, add.cell.id = "ALK08_lateE60_rep1")
ALK07_E65_rep1 <- RenameCells(object = ALK07_E65_rep1, add.cell.id = "ALK07_E65_rep1")
ALK07_lateE65_rep1 <- RenameCells(object = ALK07_lateE65_rep1, add.cell.id = "ALK07_lateE65_rep1")
ALK07_E70_rep1 <- RenameCells(object = ALK07_E70_rep1, add.cell.id = "ALK07_E70_rep1")
ALK07_earlyE75_rep1 <- RenameCells(object = ALK07_earlyE75_rep1, add.cell.id = "ALK07_earlyE75_rep1")
ALK07_lateE75_rep1 <- RenameCells(object = ALK07_lateE75_rep1, add.cell.id = "ALK07_lateE75_rep1")
ALK04_lateE75_rep2 <- RenameCells(object = ALK04_lateE75_rep2, add.cell.id = "ALK04_lateE75_rep2")
ALK07_E775_rep2 <- RenameCells(object = ALK07_E775_rep2, add.cell.id = "ALK07_E775_rep2")

# Add time stamp
ALK06_E60_rep1$"gem.group" <- "E60"
ALK06_E60_rep2$"gem.group" <- "E60"
ALK08_E60_rep3$"gem.group" <- "E60"
ALK08_lateE60_rep1$"gem.group" <- "lateE60"
ALK07_E65_rep1$"gem.group" <- "E65"
ALK07_lateE65_rep1$"gem.group" <- "lateE65"
ALK07_E70_rep1$"gem.group" <- "E70"
ALK07_earlyE75_rep1$"gem.group" <- "earlyE75"
ALK07_lateE75_rep1$"gem.group" <- "lateE75"
ALK04_lateE75_rep2$"gem.group" <- "lateE75"
ALK07_E775_rep2$"gem.group" <- "E775"

# Merge samples
v1.1.0.WT_earlyDATASETS <- merge(ALK06_E60_rep1, c(ALK06_E60_rep2,ALK08_E60_rep3,ALK08_lateE60_rep1,ALK07_E65_rep1,ALK07_lateE65_rep1,ALK07_E70_rep1,ALK07_earlyE75_rep1,ALK07_lateE75_rep1,ALK04_lateE75_rep2,ALK07_E775_rep2))
saveRDS(v1.1.0.WT_earlyDATASETS, file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS.RDS")

# Percent mito
v1.1.0.WT_earlyDATASETS[["percent.mt"]] <- PercentageFeatureSet(v1.1.0.WT_earlyDATASETS, pattern = "^mt-")

# Visualize QC metrics as a violin plot
Vln1 <- VlnPlot(v1.1.0.WT_earlyDATASETS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize QC metrics as a 2d scatter plot
plot1 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pdf(file="./0-v1.1.0/v1.1.0_QC_Plots.pdf", width=9, height=10)
print(Vln1)
print(plot1)
print(plot2)
dev.off()

# 1. Filter cells by percent.mito 
v1.1.0.WT_earlyDATASETS <- subset(v1.1.0.WT_earlyDATASETS, subset = percent.mt < 10)
plot1 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Comb1 <- CombinePlots(plots = list(plot1, plot2))

# 2. Filter cells by nCount_RNA (UMI counts)
v1.1.0.WT_earlyDATASETS <- subset(v1.1.0.WT_earlyDATASETS, subset = nCount_RNA < 50000)
plot1 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Comb2 <- CombinePlots(plots = list(plot1, plot2))


# 3. Filter cells by nFeature_RNA (measured genes)
v1.1.0.WT_earlyDATASETS <- subset(v1.1.0.WT_earlyDATASETS, subset = nFeature_RNA > 200 & nFeature_RNA < 6300)
plot1 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(v1.1.0.WT_earlyDATASETS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Comb3 <- CombinePlots(plots = list(plot1, plot2))

pdf(file="./0-v1.1.0/v1.1.0_QC_Plots-2.pdf", width=9, height=10)
print(Comb1)
print(Comb2)
print(Comb3)
dev.off()

saveRDS(v1.1.0.WT_earlyDATASETS, file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS_pre-CCregr.RDS")
#v1.1.0.WT_earlyDATASETS <- readRDS(file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS_pre-CCregr.RDS")

# cell cycle regression 
cc.genes <- read.table("./0-R_codes/2-Ressources/Cell-Cycle_genes_scRNA-Seq_MusMusculus.txt", sep="\t", header = T,fill = TRUE)
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[cc.genes$CellCyclePhase == "S",1]
g2m.genes <- cc.genes[cc.genes$CellCyclePhase == "G2/M",1]

# assign cell cycle scores
v1.1.0.WT_earlyDATASETS <- CellCycleScoring(v1.1.0.WT_earlyDATASETS, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = v1.1.0.WT_earlyDATASETS@meta.data)

##
saveRDS(v1.1.0.WT_earlyDATASETS, file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS_pre-SCT.RDS")
#v1.1.0.WT_earlyDATASETS <- readRDS(file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS_pre-SCT.RDS")

#======================================================================================================================================================================#
################################################========      Global prepared dataset analysis ## v1.2.6      ========##################################################
#======================================================================================================================================================================#

# Loading data
v1.2.6.WT_earlyDATASETS <- readRDS(file = "./0-v1.1.0/v1.1.0.WT_earlyDATASETS_pre-SCT.RDS")
v1.2.6.WT_earlyDATASETS[["Library.ident"]] <- str_split_fixed(colnames(v1.2.6.WT_earlyDATASETS), "_",2)[,1]
v1.2.6.WT_earlyDATASETS <- SCTransform(v1.2.6.WT_earlyDATASETS, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE)
v1.2.6.WT_earlyDATASETS <- RunPCA(v1.2.6.WT_earlyDATASETS, verbose = TRUE)

pdf(file="./0-v1.2.6/v1.2.6_ElbowPlot-PCA.pdf", width=9, height=10)
ElbowPlot(v1.2.6.WT_earlyDATASETS, ndims = 50, reduction = "pca")
dev.off()

print("Data loaded!")

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
v1.2.6.WT_earlyDATASETS_harmony <- RunHarmony(v1.2.6.WT_earlyDATASETS, 
				group.by.vars = c("Library.ident"), 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

print("Integration / Batch effect correction done!")

# UMAP calculation and clustering
v1.2.6.WT_earlyDATASETS_harmony <- RunUMAP(v1.2.6.WT_earlyDATASETS_harmony, dims = 1:20, reduction = "harmony", verbose = TRUE)
v1.2.6.WT_earlyDATASETS_harmony <- FindNeighbors(v1.2.6.WT_earlyDATASETS_harmony,reduction = "harmony", dims = 1:20, verbose = TRUE)
v1.2.6.WT_earlyDATASETS_harmony_res1.0 <- FindClusters(v1.2.6.WT_earlyDATASETS_harmony, verbose = TRUE, resolution = 1.0)

pdf(file="./0-v1.2.6/v1.2.6_UMAP_WT_earlyDATASETS_harmony_res1.0.pdf", width=9, height=10)
DimPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, reduction = "umap", label = T)
dev.off()

pdf(file="./0-v1.2.6/v1.2.6_UMAP_orig.ident_WT_earlyDATASETS_harmony_res1.0.pdf", width=9, height=10)
DimPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, reduction = "umap", group.by = "orig.ident")
dev.off()

pdf(file="./0-v1.2.6/v1.2.6_UMAP_Library.ident_WT_earlyDATASETS_harmony_res1.0.pdf", width=9, height=10)
DimPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, reduction = "umap", group.by = "Library.ident")
dev.off()

print("UMAP calculation and clustering done!")

# Time Stamp UMAP
pdf(file="./0-v1.2.6/v1.2.6_TimeStamp-UMAP_WT_earlyDATASETS_harmony_res1.0.pdf", width=15, height=20)
for(Day in unique(v1.2.6.WT_earlyDATASETS_harmony_res1.0@meta.data$gem.group)){
cells <- WhichCells(v1.2.6.WT_earlyDATASETS_harmony_res1.0, expression = gem.group==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, reduction = "umap", group.by = "gem.group", cells.highlight=Cells, cols.highlight=c("#4472C4"), cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

print("Time Stamp UMAP done!")

# Dot plot and Feature plots QCs
gene_list <- c("Isl1","Osr1","Tbx1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Nkx2-5","Tbx5","Gata4","Gli1","Foxf1","Hhip","Mab21l2","Hand1","Smarcd3","Shh","Ihh","Dhh")

pdf(file = "./0-v1.2.6/v1.2.6_MarkerGenes_WT_earlyDATASETS_harmony_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = i) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()

gene_list <- c("Fgf8","Isl1","Tbx1","Fgf10","Mef2c", #aSHF
			   "Aldh1a2","Nr2f2","Gata6","Hoxa1","Hoxb1","Osr1", #pSHF
			   "Gata5","Sema3d","Sfrp1","Tbx18","Tcf21","Upk3b","Wt1", #Epicardium
			   "Hand1","Mab21l2","Tdgf1","Hcn4","Hoxd1","Msx1","Myl7","Phlda2", #FHF
			   "Irx4","Wnt2","Cited1","Fbxo32","Fitm1","Gja1","Stk39","Tbx20","Tmem163", #LV
			   "Bmp4","Rspo3", #OFT 
			   "Foxc2","Foxd1","Meox1","Twist1", #PM
			   "Sema3c","Ap3b2","Arl4a","Cacna2d2","Cck","Cox6a2","Igfbpl1","Irx1","Irx2","Nppa","Pln","Ppp1r17", # RV
			   "Pou5f1","Nanog","T","Eomes","Mesp1","Pdgfra","Kdr", # Kdr = Flk1
			   "Ttn","Nkx2-5","Tbx5","Gata4",
			   "Gli1","Foxf1","Hhip",
			   "Smarcd3")

pdf(file = "./0-v1.2.6/v1.2.6_InterestingGenes_WT_earlyDATASETS_harmony_res1.0.pdf", width = 11, height = 11)
DotPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = gene_list, col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend()
dev.off()

print("Dot plot and Feature plots QCs done!")

Atlas_list <- c("Dnmt3b","Pou5f1","Epcam","Utf1",#Epiblast
			   "Nanog","Eomes", #Primitive Streak
			   "Nkx2-1", #Caudal Epiblast
			   "Ifitm3","Dnd1","Dppa3", #PGC
			   "Hhex","Otx2","Mixl1","Gsc","Lhx1","Fgf8","Pax7", #Anterior Primitive Streak
			   "Foxa2","T","Chrd","Noto", #Notochord
			   "Cer1", #Def. endoderm
			   "Apela", # Gut
			   "Mesp1","Lefty2","Mesp2", #Nascent mesoderm
			   "Osr1", #Intermediate mesoderm
			   "Cdx2","Hes7","Hoxb1","Cdx1","Gbx2", #Caudal mesoderm
			   "Tbx1","Meox1","Tcf15", #Paraxial mesoderm
			   "Aldh1a2","Dll1","Tbx6", #Somitic mesoderm 
			   "Isl1","Tcf21",#Pharyngeal mesoderm
			   "Smarcd3","Acta2","Tagln","Myl7","Nkx2-5","Myl4","Tnnt2", #Cardiomyocytes
			   "Hoxa10","Hoxa11","Tbx4", #Allantois
			   "Bmp4","Ahnak","Pmp22","Krt18","Krt8","Igf2", #Mesenchyme
			   "Etv2", #Hematoendothelial progenitors
			   "Kdr","Anxa5","Pecam1","Cdh5", #Endothelium
			   "Lmo2","Runx1", #Blood progenitors
			   "Cited4","Gata1","Hbb-bh1","Gypa","Hba-a1","Hba-a2", #Erythroid
			   "Hes3","Epha5","Cdx4","Hoxb9", #Neuro Mesodermal Progenitors
			   "Sox2","Irx3","Hesx1","Six3", #Rostral neurectoderm
			   "Sox9","Pax3","Tfap2a","Foxd3","Dlx2","Sox10", #Neural crest
			   "En1","Pax2", #Forebrain/Midbrain/Hindbrain
			   "Pax6", #Spinal cord
			   "Foxg1","Trp63","Grhl2","Grhl3", #Surface ectoderm
			   "Dkk1","Krt19","Amot","Spink1", #Visceral endoderm
			   "Emb","Cystm1","Apoe","Apoa2","Ttr", #ExE endoderm
			   "Tfap2c","Ascl2","Elf5", #ExE ectoderm
			   "Sparc","Plat","Lamb1" #Parietal endoderm
			   )
pdf(file = "./0-v1.2.6/v1.2.6_AtlasGenes_WT_earlyDATASETS_harmony_res1.0.pdf", width = 11, height = 11)
DotPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = Atlas_list, col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend()
dev.off()
 

# Violin plot QC
Vln2.1 <- VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, c("nFeature_RNA"), group.by = "orig.ident",  pt.size = 0.1)
Vln2.2 <- VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, c("nFeature_SCT"), group.by = "orig.ident",  pt.size = 0.1)
Vln3.1 <- VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, c("nCount_RNA"), group.by = "orig.ident",  pt.size = 0.1)
Vln3.2 <- VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, c("nCount_SCT"), group.by = "orig.ident",  pt.size = 0.1)
Comb5.1 <- CombinePlots(plots = list(Vln2.1, Vln2.2))
Comb5.2 <- CombinePlots(plots = list(Vln3.1, Vln3.2))

pdf(file="./0-v1.2.6/v1.2.6_VlnPlot2_PostSCT.pdf", width=9, height=10)
print(Comb5.1)
print(Vln2.1)
print(Vln2.2)
print(Comb5.2)
print(Vln3.1)
print(Vln3.2)
dev.off()

print("Vln QC done!")

# Version data save
saveRDS(v1.2.6.WT_earlyDATASETS_harmony_res1.0, file = "./0-v1.2.6/v1.2.6.WT_earlyDATASETS_harmony_res1.0.RDS")
#v1.2.6.WT_earlyDATASETS_harmony_res1.0 <- readRDS(file = "./0-v1.2.6/v1.2.6.WT_earlyDATASETS_harmony_res1.0.RDS")

#======================================================================================================================================================================#
########################################========      Mesp1 Lineage & Expressing dataset construction ## v1.2.6.5      ========#########################################
#======================================================================================================================================================================#
v1.2.6.WT_earlyDATASETS_harmony_res1.0 <- readRDS(file = "./0-v1.2.6/v1.2.6.WT_earlyDATASETS_harmony_res1.0.RDS")

#Ai14-transgene = Mesp1-cre lineage / eGFP-transgene = Smarcd3+ cardiac precursor cells

pdf(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/Ai14-transgene_eGFP-transgene_v1.2.6.WT_earlyDATASETS_harmony_res1.0.pdf", width = 10, height = 10)
print(FeaturePlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = 'Ai14-transgene') & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
print(VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, 'Ai14-transgene', pt.size = 0.1))
print(FeaturePlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = 'eGFP-transgene') & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
print(VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, 'eGFP-transgene', pt.size = 0.1))
dev.off()

pdf(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/Mesp1_v1.2.6.WT_earlyDATASETS_harmony_res1.0.pdf", width = 10, height = 10)
print(FeaturePlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, features = 'Mesp1') & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
print(VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, 'Mesp1', pt.size = 0.1))
print(VlnPlot(v1.2.6.WT_earlyDATASETS_harmony_res1.0, 'Mesp1', pt.size = 0.1)+ stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95))
dev.off()

###################################################################
##  Subset 1  # clusters to keep: 0,4,9,12,14,15,17,20,21,24,33  ##
###################################################################
v1.2.6.5_SUBSET1 <- subset(v1.2.6.WT_earlyDATASETS_harmony_res1.0, idents = c(0,4,9,12,14,15,17,20,21,24,33), invert=FALSE)
v1.2.6.5_SUBSET1 <- SCTransform(v1.2.6.5_SUBSET1, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE)
v1.2.6.5_SUBSET1 <- RunPCA(v1.2.6.5_SUBSET1, verbose = TRUE)

pdf(file="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_ElbowPlot-PCA.pdf", width=9, height=10)
ElbowPlot(v1.2.6.5_SUBSET1, ndims = 50, reduction = "pca")
dev.off()

print("Data loaded!")

# 21,643 cells
# 30 PCs

# Harmony integration / Batch effect correction
#Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity
v1.2.6.5_SUBSET1_harmony <- RunHarmony(v1.2.6.5_SUBSET1, 
				group.by.vars = c("Library.ident"), 
				reduction = "pca", 
				assay.use = "SCT", 
				reduction.save = "harmony")

print("Integration / Batch effect correction done!")

# UMAP calculation and clustering
v1.2.6.5_SUBSET1_harmony <- RunUMAP(v1.2.6.5_SUBSET1_harmony, dims = 1:30, reduction = "harmony", verbose = TRUE)
v1.2.6.5_SUBSET1_harmony <- FindNeighbors(v1.2.6.5_SUBSET1_harmony,reduction = "harmony", dims = 1:30, verbose = TRUE)
v1.2.6.5_SUBSET1_harmony_res1.0 <- FindClusters(v1.2.6.5_SUBSET1_harmony, verbose = TRUE, resolution = 1.0)

pdf(file="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/UMAP_v1.2.6.5_SUBSET1_harmony_res1.0.pdf", width=9, height=10)
print(DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", label = T))
print(DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", group.by = "orig.ident"))
print(DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", group.by = "Library.ident"))
dev.off()

print("UMAP calculation and clustering done!")

# Time Stamp UMAP
pdf(file="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/TimeStamp-UMAP_v1.2.6.5_SUBSET1_harmony_res1.0.pdf", width=15, height=20)
for(Day in unique(v1.2.6.5_SUBSET1_harmony_res1.0@meta.data$gem.group)){
cells <- WhichCells(v1.2.6.5_SUBSET1_harmony_res1.0, expression = gem.group==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", group.by = "gem.group", cells.highlight=Cells, cols.highlight=c("#4472C4"), cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

print("Time Stamp UMAP done!")

# Dot plot and Feature plots QCs
gene_list <- c("Fgf5","Isl1","Osr1","Tbx1","Myl7","Ttn","Pou5f1","Nanog","T","Eomes","Mesp1","Asb2","Nkx2-5","Tbx5","Gata4","Gata5","Gata6","Gli1","Foxf1","Hhip","Ptch1","Fgf8","Fgf10","Mab21l2","Hand1","Smarcd3","Ai14-transgene","eGFP-transgene","Myl2","Hcn4","Hopx")
pdf(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/MarkerGenes_v1.2.6.5_SUBSET1_harmony_res1.0.pdf", width = 10, height = 10)
for(i in gene_list){			   
print(FeaturePlot(v1.2.6.5_SUBSET1_harmony_res1.0, features = i) & scale_color_gradientn(colors = plasma(n = 10, direction = -1)))
}	   
dev.off()
# Version data save
saveRDS(v1.2.6.5_SUBSET1_harmony_res1.0, file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")
#v1.2.6.5_SUBSET1_harmony_res1.0 <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")

print("Version saved!")

################################
## Top markers identification ##
################################
#v1.2.6.5_SUBSET1_harmony_res1.0 <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")

# Find markers
v1.2.6.5_SUBSET1_harmony_res1.0 <- PrepSCTFindMarkers(v1.2.6.5_SUBSET1_harmony_res1.0)
v1.2.6.5_SUBSET1_harmony_res1.0.markers <- FindAllMarkers(v1.2.6.5_SUBSET1_harmony_res1.0, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.6,recorrect_umi = FALSE)
write.csv(v1.2.6.5_SUBSET1_harmony_res1.0.markers, file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/TopMarkers/v1.2.6.5_SUBSET1_harmony_res1.0_all-markers.csv")

# Marker heatmap
v1.2.6.5_SUBSET1_harmony_res1.0.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15
	
pdf(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/TopMarkers/TopMarkers_v1.2.6.5_SUBSET1_harmony_res1.0.pdf", width = 11, height = 11)
print(DoHeatmap(v1.2.6.5_SUBSET1_harmony_res1.0, features = top15$gene) + NoLegend())
dev.off()

pdf(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/TopMarkers/DotPlot_v1.2.6.5_SUBSET1_harmony_res1.0.pdf", width = 15, height = 40)
print(DotPlot(v1.2.6.5_SUBSET1_harmony_res1.0, features = unique(top15$gene), col.max = 3, col.min = -2, cols = c("grey","red"), assay = "SCT", dot.scale = 6) +coord_flip() + RotatedAxis() + xlab('') + ylab('') + NoLegend())
dev.off()

##########################################
## SingleR - Mouse gastrulation dataset ##
##########################################
# sbatch script used
library(SingleR)

v1.2.6.5_SUBSET1_harmony_res1.0 <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")

ref <- readRDS("/project/imoskowitz/wiktor/4-singlecell_reference_datasets/RefEmbryoAtlasData.rds")
ref <- ref[, ref$doublet == FALSE] # remove doublets
ref <- ref[, !is.na(ref$celltype)] # remove empty cell types
ref_stage <- ref[, ref$stage %in% c("E6.5","E6.75","E7.0","E7.25","E7.5","E7.75")] # keep only stage e8.25 - Keep the time stage you want
ref_stage <- scuttle::logNormCounts(ref_stage) # singleR needs log normalized reference
sce_v1.2.6.5.SUBSET1 <- as.SingleCellExperiment(v1.2.6.5_SUBSET1_harmony_res1.0)

# Convert symbol to ensembl
require(org.Mm.eg.db)
ens <- mapIds(org.Mm.eg.db,
	        keys = rownames(sce_v1.2.6.5.SUBSET1),
		  column = 'ENSEMBL',
		  keytype = 'SYMBOL')
all(rownames(sce_v1.2.6.5.SUBSET1) == names(ens))

keep <- !is.na(ens)
ens <- ens[keep]
sce_v1.2.6.5.SUBSET1 <- sce_v1.2.6.5.SUBSET1[keep,]
rownames(sce_v1.2.6.5.SUBSET1) <- ens

# Run single R
mouse_gastrulation_label <- SingleR(test=sce_v1.2.6.5.SUBSET1, ref = ref_stage, labels = ref_stage$celltype, clusters = sce_v1.2.6.5.SUBSET1$seurat_clusters)
saveRDS(mouse_gastrulation_label, file ="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/SingleR_MouseGastrulationAtlas/mouse_gastrulation_label.RDS")

v1.2.6.5_SUBSET1_harmony_res1.0$mouse_gastrulation_label <- factor(v1.2.6.5_SUBSET1_harmony_res1.0$seurat_clusters, labels = mouse_gastrulation_label$labels)
saveRDS(v1.2.6.5_SUBSET1_harmony_res1.0, file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/SingleR_MouseGastrulationAtlas/v1.2.6.5_SUBSET1_harmony_res1.0_SlingleR.RDS")
#v1.2.6.5_SUBSET1_harmony_res1.0_SlingleR <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/SingleR_MouseGastrulationAtlas/v1.2.6.5_SUBSET1_harmony_res1.0_SlingleR.RDS")

pdf(file="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/SingleR_MouseGastrulationAtlas/UMAP_v1.2.6.5_SUBSET1_harmony_res1.0_SlingleR.pdf", width=9, height=10)
print(DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", group.by = "mouse_gastrulation_label"))
print(DimPlot(v1.2.6.5_SUBSET1_harmony_res1.0, reduction = "umap", group.by = "mouse_gastrulation_label",label = T)+NoLegend())
dev.off()

pdf(file="./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/SingleR_MouseGastrulationAtlas/UMAP_v1.2.6.5_SUBSET1_harmony_res1.0_SlingleR-plotScoreHeatmap.pdf", width=9, height=10)
print(plotScoreHeatmap(mouse_gastrulation_label, clusters=rownames(mouse_gastrulation_label)))
dev.off()

#======================================================================================================================================================================#
########################################========            Subset 1 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
v1.2.6.5_SUBSET1_harmony_res1.0 <- readRDS(file = "./Krup_GSE208153_scRNAseq/0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0.RDS")
CellsToRemove_Krup <- read.csv("./Krup-Tyser-Pijuan_Joint/1-Merge1/PrepSub1_CellsToRemove_Krup.csv")

v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan <- subset(v1.2.6.5_SUBSET1_harmony_res1.0, cells = CellsToRemove_Krup$x, invert=TRUE)
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan <- SCTransform(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan, file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan.RDS")
#v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan.RDS")

#======================================================================================================================================================================#
########################################========            Subset 2 after merged atlas analysis           ========#########################################
#======================================================================================================================================================================#
	# In addition to previous 'subset1'
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan.RDS")
CellsToRemove2_Krup <- read.csv("./Krup-Tyser-Pijuan_Joint/2-Merge2/PrepSub2_CellsToRemove_Krup.csv")

v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan <- subset(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge_Krup_Tyser_Pijuan, cells = CellsToRemove2_Krup$x, invert=TRUE)
v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan <- SCTransform(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "nCount_RNA","nFeature_RNA"), verbose = TRUE)

# Version data save
saveRDS(v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan, file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan.RDS")
#v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan <- readRDS(file = "./0-v1.2.6/v1.2.6.5_Mesp1Lineage_Mesp1Expressing/1-Subset1/v1.2.6.5_SUBSET1_harmony_res1.0_subForMerge3_Krup_Tyser_Pijuan.RDS")

#2254 features across 11975 samples within 2 assays
#Active assay: SCT (19413 features, 3000 variable features)
