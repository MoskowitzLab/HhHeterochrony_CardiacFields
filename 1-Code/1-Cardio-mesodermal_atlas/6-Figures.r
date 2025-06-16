### Figure code v4
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

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- PrepSCTFindMarkers(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0)

## Color codes ##
AnnotationColors = c("#abdda4","#f46d43","#4178ae","#af5b39","#2095de","#fb9ac2","#E41A1C","#a9a9a9","#ffce00","#40928a","#5e4fa2")
names(AnnotationColors) = c("Nascent mesoderm","Mixed mesoderm","Primitive streak","Mesenchyme","ExE mesoderm","JCF/Mesenchyme","CMs","Paraxial/Pharyngeal mesoderm","Differentiating CMs","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm")

TimepointColors = c("#F58484","#DB7297","#C261A5","#A451A8","#70428F","#443475","#272C5C","#1B2C42","#102429","#060F0E")
names(TimepointColors) = c("E6.0","E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")

##############################################################################################################################################
############################################                       Figure 1                       ############################################
##############################################################################################################################################

### Annotated UMAP (Res1) 
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters %in% c("0","15","17") ~ "Primitive streak",
	seurat_clusters %in% c("3","10","16") ~ "Nascent mesoderm",
	seurat_clusters %in% c("1","14") ~ "Mixed mesoderm",
	seurat_clusters %in% "13" ~ "Paraxial/Pharyngeal mesoderm",
	seurat_clusters %in% c("5","8") ~ "ExE mesoderm",
	seurat_clusters %in% "6" ~ "aSHF/Pharyngeal mesoderm",
	seurat_clusters %in% "9" ~ "pSHF/Pharyngeal mesoderm",
	seurat_clusters %in% c("2","11") ~ "Mesenchyme",
	seurat_clusters %in% "4" ~ "JCF/Mesenchyme",
	seurat_clusters %in% "12" ~ "Differentiating CMs",
	seurat_clusters %in% "7" ~ "CMs",
	TRUE ~ seurat_clusters
  ))
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation <- factor(
  Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation,
  levels = c("Primitive streak", "Nascent mesoderm","Mixed mesoderm","Paraxial/Pharyngeal mesoderm","ExE mesoderm","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm","Mesenchyme","JCF/Mesenchyme","Differentiating CMs","CMs")
)

pdf(file="./3-Merge3/Figures/Figure1_UMAP_Annotation.pdf", width=10, height=10)
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T, pt.size = 0.7, group.by = "Annotation", raster = FALSE, cols = AnnotationColors))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = F, pt.size = 0.7, group.by = "Annotation", raster = FALSE, cols = AnnotationColors))
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = T, pt.size = 0.7, group.by = "Annotation", raster = FALSE, cols = AnnotationColors) + NoLegend())
print(DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", label = F, pt.size = 0.7, group.by = "Annotation", raster = FALSE, cols = AnnotationColors) + NoLegend())
dev.off()

### Heatmap annotated UMAP
AnnoGenes <- c(
			"Nanog","Eomes","T",# Primitive streak
			"Mesp1","Mesp2","Lefty2","Tdgf1",# Nascent mesoderm
			"Evx1",# Mixed mesoderm
			"Tcf15",# Paraxial/Pharyngeal mesoderm
			"Cdx1","Cdx2","Cdx4","Bmp4",# ExE mesoderm
			"Isl1","Tbx1","Fgf10",# aSHF/Pharyngeal mesoderm
			"Osr1","Aldh1a2","Arg1",# pSHF/Pharyngeal mesoderm
			"Igf2","Krt8","Krt18","Pmp22","Ahnak",# Mesenchyme
			"Hand1","Foxf1","Mab21l2",# JCF/Mesenchyme
			"Ttn","Tnnt2","Myl7","Acta2"# Cardiomyocytes
)

Idents(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0) <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation

expression_matrix <- GetAssayData(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, assay = "SCT", slot = "data")
annotation_labels <- setNames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation, rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data))

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
  col = list(Annotation = AnnotationColors),
  annotation_label = "Cell type",
  show_legend = FALSE 
)

pdf(file = "./3-Merge3/Figures/Figure1_Heatmap_AveragedExpression.pdf", width = 10, height = 6.9)
Heatmap(
  matrix = as.matrix(scaled_expression),
  name = NULL,
  width = unit(4.6, "in"),
  col = colorRamp2(c(-2,0,3), hcl_palette = "viridis"),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = NULL,
  row_title = NULL,
  column_names_rot = 45,
  bottom_annotation = bottom_annotation,
  heatmap_legend_param = list(
    title = "Expression levels",
    at = c(-2, 3),
    labels = c("Low", "High"),
    title_position = "topcenter",
    labels_gp = gpar(fontsize = 8),
    legend_direction = "vertical",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    border = FALSE
  )
)
dev.off()

### Ribbon plot - Contribution of different cell type to stages
library(reshape)

meta <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data
meta$Annotation <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation
timepointOrder <- c("E6.0","E6.5","E6.75","E7.0","E7.25","E7.5","E7.75","E8.0","E8.25","E8.5")

contrib_known2 = table(meta$TimeStampMerged, meta$Annotation)/rowSums(table(meta$TimeStampMerged, meta$Annotation))

pdf2 = melt(contrib_known2)
names(pdf2) = c("stage", "cluster", "fraction")

pdf2$cluster <- factor(pdf2$cluster, levels = levels(meta$Annotation))
AnnotationColors <- AnnotationColors[levels(meta$Annotation)] 

max_positions <- pdf2 %>%
filter(!stage %in% c("E6.0", "E8.5")) %>% #avoid selecting borders for plotting
  group_by(cluster) %>%
  slice(which.max(fraction))

pdf2$cluster_numeric <- as.numeric(as.character(pdf2$cluster))
max_positions$cluster_numeric <- as.numeric(as.character(max_positions$cluster))

Data2 <- NULL
for(clus in unique(max_positions$cluster)){
temp <- pdf2[pdf2$stage == max_positions$stage[max_positions$cluster == clus],]
temp2 <- temp[!(temp$cluster_numeric <= as.numeric(clus)),]
sumFrac <- sum(temp2$fraction)
yVal <- sumFrac + (temp$fraction[temp$cluster == clus])/2
Data <- data.frame(cluster = clus, y_coord = yVal, stage = max_positions$stage[max_positions$cluster == clus])
Data2 <- rbind(Data2,Data)
}

pdf(file=paste0("./3-Merge3/Figures/Figure1_Ribbon_ContribTypetoTime.pdf"), width=7, height=7)
print(ggplot(pdf2, aes(x = factor(stage, ordered = TRUE, levels = timepointOrder), y = fraction, fill = cluster, group = cluster)) +
  geom_area(position = "stack") +
  geom_text(data = Data2, aes(x = factor(stage, ordered = TRUE, levels = timepointOrder), y = y_coord, label = cluster), size = 6, color = "black") +
  scale_fill_manual(values = AnnotationColors, labels = names(AnnotationColors), name = "") +
  labs(x = NULL, y= "Cell fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.background = element_blank(),panel.grid.minor.y = element_blank(),axis.line.y = element_line(color = "black"))+
  scale_y_continuous(expand = c(0, 0))
)
print(ggplot(pdf2, aes(x = factor(stage, ordered = TRUE, levels = timepointOrder), y = fraction, fill = factor(cluster, levels = sort(unique(meta$Annotation))), group = cluster)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = AnnotationColors, labels = names(AnnotationColors), name = "") +
  labs(x = NULL, y= "Cell fraction") +
  theme(aspect.ratio = 1,axis.text.x = element_text(size = 14, angle = 45, hjust = 1),axis.text.y = element_text(size = 14),panel.background = element_blank(),panel.grid.minor.y = element_blank(),axis.line.y = element_line(color = "black"))+
  scale_y_continuous(expand = c(0, 0)) +
  NoLegend()
)
print(ggplot(pdf2, aes(x = factor(stage, ordered = TRUE, levels = timepointOrder), y = fraction, fill = factor(cluster, levels = sort(unique(meta$Annotation))), group = cluster)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = AnnotationColors, labels = names(AnnotationColors), name = "") +
  labs(x = NULL, y= "Cell fraction") +
  theme(aspect.ratio = 1,axis.text.x = element_text(size = 14, angle = 45, hjust = 1),axis.text.y = element_text(size = 14),panel.background = element_blank(),panel.grid.minor.y = element_blank(),axis.line.y = element_line(color = "black"))+
  scale_y_continuous(expand = c(0, 0)) +
  NoLegend() +
  annotate("text", x = 3.5, y = 0.88, label = "Primitive streak", size = 5, color = "white") +
  annotate("text", x = 3, y = 0.58, label = "Nascent mesoderm", size = 5, color = "white") +
  annotate("text", x = 3, y = 0.25, label = "Mixed mesoderm", size = 5, color = "white") +
  annotate("text", x = 7.5, y = 0.89, label = "Paraxial mesoderm", size = 5, color = "white") +
  annotate("text", x = 8.5, y = 0.82, label = "ExE mesoderm", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.68, label = "aSHF", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.55, label = "pSHF", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.4, label = "Mesenchyme", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.26, label = "JCF", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.18, label = "Diff. CMs", size = 5, color = "white") +
  annotate("text", x = 9, y = 0.08, label = "CMs", size = 5, color = "white")
)
dev.off()

### Timepoint UMAP
pdf(file="./3-Merge3/Figures/Figure1-TimeStampMerged.pdf", width=10, height=10)
for(Day in unique(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimeStampMerged)){
cells <- WhichCells(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimeStampMerged ==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", group.by = "TimeStampMerged", cells.highlight=Cells, cols.highlight=c("#4472C4"),sizes.highlight = 0.6, pt.size = 0.4, cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()

### Timepoint UMAP (compacted version)
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(TimePointCompacted = case_when(
    TimeStampMerged %in% "E6.0" ~ "Early",
	TimeStampMerged %in% "E6.5" ~ "Early",
	TimeStampMerged %in% "E6.75" ~ "Early",
	TimeStampMerged %in% "E7.0" ~ "Early",
	TimeStampMerged %in% "E7.25" ~ "Early",
	TimeStampMerged %in% "E7.5" ~ "Switch",
	TimeStampMerged %in% "E7.75" ~ "Late",
	TimeStampMerged %in% "E8.0" ~ "Late",
	TimeStampMerged %in% "E8.25" ~ "Late",
	TimeStampMerged %in% "E8.5" ~ "Late",
	TRUE ~ TimeStampMerged
  ))

pdf(file="./3-Merge3/Figures/Figure1-TimeStampMerged-Compacted.pdf", width=10, height=10)
for(Day in unique(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimePointCompacted)){
cells <- WhichCells(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimePointCompacted ==Day)
Cells <- setNames(list(cells), Day)
c <- DimPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, reduction = "umap", group.by = "TimePointCompacted", cells.highlight=Cells, cols.highlight=c("#4472C4"),sizes.highlight = 0.6, pt.size = 0.4, cols=c("#E4E4E4"))+
  ggtitle(Day)+
  NoLegend()
print(c)
}
dev.off()


### Hh score Vlnplot per cell type
UCell_score <- read.csv("./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0[["Limited_Hh_list_UCell"]] <- UCell_score$UCell[match(rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), UCell_score$X)]
Idents(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0) <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Annotation

median_data <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  group_by(Annotation) %>%
  summarise(median_score = median(Limited_Hh_list_UCell))

pdf(file = "./3-Merge3/Figures/Figure1_ViolinPlot_HhScore-perCellType.pdf", width = 11, height = 11)	   
print(VlnPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, "Limited_Hh_list_UCell", pt.size = 0)+ 
	    NoLegend() +
		geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
		scale_fill_manual(values = AnnotationColors) +
		labs(title = NULL,
			     x = NULL,              
                 y = "Hedgehog Score (AU)") + 
		theme(axis.title.x = element_blank()))
print(VlnPlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, "Limited_Hh_list_UCell", pt.size = 0)+ 
	    NoLegend() +
		geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
		scale_fill_manual(values = AnnotationColors) +
		labs(title = NULL,
			     x = NULL,              
                 y = "Hedgehog Score (AU)") + 
		theme(axis.title.x = element_blank(),
		axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)))
dev.off()

# Calculate percentage of cells with Limited_Hh_list_UCell > 0 per Annotation group
percentage_df <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  group_by(Annotation) %>%
  summarise(percentage = sum(Limited_Hh_list_UCell > 0) / n() * 100)
percentage_df

### Hh score UMAP
pdf(file = "./3-Merge3/Figures/Figure1_HhScore-UMAP.pdf", width = 10, height = 10)	   
print(FeaturePlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, features = "Limited_Hh_list_UCell", order = TRUE, pt.size = 0.6) + 
	  scale_color_gradientn(limits = c(0, 0.2),colors = c("#E4E4E4","#C71717"), na.value = "#F6F0F0")
)
dev.off()	

### Hh score UMAP per timepoint
timepoints <- unique(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimeStampMerged)

pdf(file = "./3-Merge3/Figures/Figure1_HhScore_ScorePerTimePoint.pdf", width = 10, height = 10)
  for (Day in timepoints) {
  cells <- WhichCells(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimeStampMerged == Day)
  
    # Create a vector with score values, setting NA for cells not in the timepoint
    expr <- FetchData(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, vars = "Limited_Hh_list_UCell")
    expr[!rownames(expr) %in% cells, "Limited_Hh_list_UCell"] <- NA
    
    # Add expression data to metadata
    Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Score_expr <- expr[, "Limited_Hh_list_UCell"]
    
    p <- FeaturePlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, features = "Score_expr", reduction = "umap", pt.size = 1.1, order = TRUE) +
      scale_color_gradientn(limits = c(0,0.2),colors = c("#A8C9D8","#f64f59","#C71717"), na.value = "#E4E4E4") +
      ggtitle(paste(Day, "-", "Limited_Hh_list_UCell"))
    print(p)
  }
dev.off()

### Hh score UMAP per timepoint (compacted version)
Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data %>%
  mutate(TimePointCompacted = case_when(
    TimeStampMerged %in% "E6.0" ~ "Early",
	TimeStampMerged %in% "E6.5" ~ "Early",
	TimeStampMerged %in% "E6.75" ~ "Early",
	TimeStampMerged %in% "E7.0" ~ "Early",
	TimeStampMerged %in% "E7.25" ~ "Early",
	TimeStampMerged %in% "E7.5" ~ "Switch",
	TimeStampMerged %in% "E7.75" ~ "Late",
	TimeStampMerged %in% "E8.0" ~ "Late",
	TimeStampMerged %in% "E8.25" ~ "Late",
	TimeStampMerged %in% "E8.5" ~ "Late",
	TRUE ~ TimeStampMerged
  ))

timepoints <- unique(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$TimePointCompacted)

pdf(file = "./3-Merge3/Figures/Figure1_HhScore_ScorePerTimePoint_Compacted.pdf", width = 10, height = 10)
  for (Day in timepoints) {
  cells <- WhichCells(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, expression = TimePointCompacted == Day)
  
    # Create a vector with score values, setting NA for cells not in the timepoint
    expr <- FetchData(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, vars = "Limited_Hh_list_UCell")
    expr[!rownames(expr) %in% cells, "Limited_Hh_list_UCell"] <- NA
    
    # Add expression data to metadata
    Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data$Score_expr <- expr[, "Limited_Hh_list_UCell"]
    
    p <- FeaturePlot(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, features = "Score_expr", reduction = "umap", pt.size = 1.1, order = TRUE) +
      scale_color_gradientn(limits = c(0,0.2),colors = c("#A8C9D8","#f64f59","#C71717"), na.value = "#E4E4E4") +
      ggtitle(paste(Day, "-", "Limited_Hh_list_UCell"))
    print(p)
  }
dev.off()


##############################################################################################################################################
############################################                       Figure 2                       ############################################
##############################################################################################################################################

#########################################
### FHF Slingshot trajectory analysis ###
#########################################
FHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/FHF2/FHF_Slingshot_SelectedCells.RDS")

FHF_Slingshot_SelectedCells@meta.data <- FHF_Slingshot_SelectedCells@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters %in% c("0","15","17") ~ "Primitive streak",
	seurat_clusters %in% c("3","10","16") ~ "Nascent mesoderm",
	seurat_clusters %in% c("1","14") ~ "Mixed mesoderm",
	seurat_clusters %in% "13" ~ "Paraxial/Pharyngeal mesoderm",
	seurat_clusters %in% c("5","8") ~ "ExE mesoderm",
	seurat_clusters %in% "6" ~ "aSHF/Pharyngeal mesoderm",
	seurat_clusters %in% "9" ~ "pSHF/Pharyngeal mesoderm",
	seurat_clusters %in% c("2","11") ~ "Mesenchyme",
	seurat_clusters %in% "4" ~ "JCF/Mesenchyme",
	seurat_clusters %in% "12" ~ "Differentiating CMs",
	seurat_clusters %in% "7" ~ "CMs",
	TRUE ~ seurat_clusters
  ))
FHF_Slingshot_SelectedCells@meta.data$Annotation <- factor(
  FHF_Slingshot_SelectedCells@meta.data$Annotation,
  levels = c("Primitive streak", "Nascent mesoderm","Mixed mesoderm","Paraxial/Pharyngeal mesoderm","ExE mesoderm","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm","Mesenchyme","JCF/Mesenchyme","Differentiating CMs","CMs")
)

# UMAP selected cells
FHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"), invert=FALSE)
FHF2_dataSub <- subset(x = FHF2_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,7,8,11,12), invert=FALSE)

good <- data.frame(cells = rownames(FHF2_dataSub@meta.data), annot = FHF2_dataSub@meta.data$Annotation)
all <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), annot = "Other")
rest <- all[!(all$cells %in% good$cells),]
full <- rbind(good,rest)

temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0

temp@meta.data$CustomColors <- full$annot[match(rownames(temp@meta.data), full$cells)]
AnnotColors <- c(AnnotationColors, "Other" = "#E4E4E4")

pdf(file="./3-Merge3/Figures/Figure2_UMAP_selectedCells_ForFHFtrajectory.pdf", width=10, height=10)
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
dev.off()

# UMAP plot pseudotime
pdf(file = "./3-Merge3/Figures/Figure2_TrajectoryViz_FHF.pdf", width = 10, height = 10)
  good2 <- data.frame(cells = colnames(FHF_Slingshot_SelectedCells), Pseudotime = FHF_Slingshot_SelectedCells@meta.data$SlingshotPseudotime)
  all2 <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), Pseudotime = NA)
  rest2 <- all2[!(all2$cells %in% good2$cells),]
  full2 <- rbind(good2,rest2)
temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
temp@meta.data$Pseudo <- full2$Pseudotime[match(rownames(temp@meta.data), full2$cells)]
print(FeaturePlot(temp, features = "Pseudo", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
dev.off()

#pdf(file="./3-Merge3/Figures/Figure2_TrajectoryViz_FHF.pdf", width=10, height=10)
#print(DimPlot(FHF_Slingshot_SelectedCells, reduction = "umap", group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(FHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 2), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(FHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(FHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(2, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(FeaturePlot(FHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "umap", order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#print(FeaturePlot(FHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "harmony",dims = c(1, 2), order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#dev.off()

# Ridgeplot timestamp contribution to pseudotime
preTable <- FHF_Slingshot_SelectedCells@meta.data[,c("TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(preTable) <- c("TimeStampMerged","Pseudotime","Annotation")
pdf(file = "./3-Merge3/Figures/Figure2_FHF_RidgePlot.pdf", width = 8, height = 8)
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 1.5, size = 0.3) +
    scale_fill_gradientn(
      colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
      name = "Pseudotime"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove grid lines
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),  # Add space for arrow and title
      axis.text.y = element_text(size = 20, margin = margin(l = 10)),  # Push timepoint labels to the right
      axis.text.x = element_text(size = 20),  # Adjust size of x-axis labels
      axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1)  # Adjust x-axis title
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)),  # First and last breaks
      labels = c("Early", "Late")  # Replace numerical x-axis labels
    ) +
    labs(x = "Pseudotime", y = "Cell density") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    ) 
 )
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges(aes(fill = Annotation), scale = 1.5, size = 0.3) +
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
# Calculate densities and interpolate to match original data length
density_data <- preTable %>%
  group_by(TimeStampMerged) %>%
  mutate(
    Density = approxfun(density(Pseudotime, adjust = 1)$x, 
                        density(Pseudotime, adjust = 1)$y)(Pseudotime)
  )

# Plot with accurate density and coloring by Annotation
print(ggplot(density_data, aes(x = Pseudotime, y = TimeStampMerged, height = Density, fill = Annotation)) +
    geom_density_ridges(stat = "identity", scale = 1, size = 0.3) + # Use stat = "identity"
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
dev.off()

# Heatmap + GO terms top 500 DEGs along the pseudotime 
FHF_DEGlist <- read.csv("./3-Merge3/Slingshot/FHF2/tradeSeq/FHF_DEGlist.csv")
rownames(FHF_DEGlist) <- FHF_DEGlist$X
FHF_DEGlist <- FHF_DEGlist[order(-FHF_DEGlist$waldStat), ]
top500_FHF <- rownames(FHF_DEGlist[1:500,])

magict4_data <- readRDS(file = "./3-Merge3/DATANormcountsAfterMAGIC_t4.RDS")
meta <- FHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
magict4_data <- magict4_data[, rownames(meta_sort)]
magict4_data <- magict4_data[top500_FHF, ]
tmagict4_data <- t(scale(t(magict4_data)))  # Row scaling (genes)

	#kmeans
set.seed(123)
wss <- sapply(1:15, 
              function(k){kmeans(tmagict4_data, k, nstart=25,iter.max = 15)$tot.withinss})
pdf(file = "./3-Merge3/Slingshot/FHF2/tradeSeq/FHF_HeatmapTop500_kmean_ElbowPlot.pdf", width = 8, height = 8)	
plot(1:15, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()
	
km_clusters <- data.frame(cluster = kmeans(tmagict4_data, centers = 6, iter.max = 15, nstart = 25)$cluster) # 25 random sets tested, 15 iteration for convergence
km_clusters$geneid <- rownames(km_clusters)
write.csv(as.data.frame(km_clusters), file = "./3-Merge3/Slingshot/FHF2/tradeSeq/Top500DEGs_km_clusters_FHF.csv")
#km_clusters <- read.csv("./3-Merge3/Slingshot/FHF2/tradeSeq/Top500DEGs_km_clusters_FHF.csv", row.names = 1)

row_order = hclust(dist(tmagict4_data))

	#GO terms
universe <- rownames(FHF_Slingshot_SelectedCells)
FHF_go <- compareCluster(geneid~cluster, data = km_clusters, fun = 'enrichGO', OrgDb  = "org.Mm.eg.db", ont='BP', qvalueCutoff = 0.05,  keyType = "SYMBOL", universe = universe) #readable = T,
write.csv(as.data.frame(FHF_go@compareClusterResult), file = "./3-Merge3/Slingshot/FHF2/tradeSeq/Top500DEGs_GOterms_FHF.csv")
#FHF_go <- read.csv("./3-Merge3/Slingshot/FHF2/tradeSeq/Top500DEGs_GOterms_FHF.csv", row.names = 1)

		# extract top 5 terms per cluster
n <- 5
FHF_go_list <- FHF_go@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(qvalue, n = n) %>%  # Select top N based on qvalue
  summarise(GO_Terms = list(Description), .groups = "drop") %>%
  deframe()

	#Heatmap
		# Annotate go terms on the heatmap
		# Create annotation object using anno_textbox()
go_anno <- rowAnnotation(
  textbox = anno_textbox(
    km_clusters$cluster, 
    FHF_go_list, 
    #by = 'anno_block',
    word_wrap = T, 
    add_new_line = T,
   background_gp = gpar(fill = 'transparent')
  ),
  width = unit(4.5, "in")
)
combined_annotation <- HeatmapAnnotation(
  pseudotime = meta_sort$dptval,
  Annotation = meta_sort$Annotation,
  TimeStampMerged = meta_sort$TimeStampMerged,
  col = list(
    pseudotime = colorRamp2(
      seq(min(meta_sort$dptval), max(meta_sort$dptval), length.out = 256),  # Divide pseudotime into 256 gradient steps
      viridis(256, option = "plasma")  # Full plasma palette with 256 colors
    ),
    Annotation = AnnotationColors,
    TimeStampMerged = TimepointColors
  ),
  annotation_legend_param = list(
    pseudotime = list(title = "Pseudotime", at = c(min(meta_sort$dptval), max(meta_sort$dptval)),labels = c("Early", "Late")),
    Annotation = list(title = "Cell Type"),
    TimeStampMerged = list(title = "Timepoints")
  ),
  annotation_label = list(pseudotime = "Pseudotime", Annotation = "Cell type", TimeStampMerged = "Timepoints")
)
pdf(file = "./3-Merge3/Figures/Figure2_FHF_HeatmapTop500.pdf", width = 26, height = 8)	   
draw(Heatmap(
  tmagict4_data,
  cluster_columns = FALSE,  # No clustering for cells
  cluster_rows = TRUE,      # Cluster genes
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE, 
  show_column_names = FALSE,
  col = colorRamp2(c(-1.5,0,1.5), hcl_palette = "viridis"),
  row_split = km_clusters$cluster,
  row_gap = unit(1, "mm"),
  row_title = NULL,
  right_annotation = go_anno,
  heatmap_legend_param = list(
    title = "Expression levels",
    at = c(-2, 2),
    labels = c("Low", "High")),
    width = unit(16.6, "in")
),
annotation_legend_side = "right",
heatmap_legend_side = "right"
)
dev.off()

## Plot unique gene along Pseudotime ## Curve version2 =========
	# logCPM
rna_data <- FHF_Slingshot_SelectedCells[["RNA"]]@data
meta <- FHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
rna_data <- as.matrix(rna_data[, rownames(meta_sort)])

GENE_LIST <- c("Mesp1","Gata4","Nkx2-5","Tbx5","Myl7","Asb2","Ttn","Tnnt2")

pdf(file = "./3-Merge3/Figures/Figure2_FHF_GeneCurve.pdf", width = 8, height = 8)
for(GENE in GENE_LIST){
ConvertTable <- NULL
ConvertTable <- data.frame(logCPM = rna_data[rownames(rna_data) == GENE,], row.names = rownames(t(rna_data)))
ConvertTable$dptval <- meta_sort$dptval[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$Annotation <- meta_sort$Annotation[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$TimeStampMerged <- meta_sort$TimeStampMerged[match(rownames(ConvertTable), rownames(meta_sort))]
    print(ggplot(ConvertTable, aes(x = dptval, y = logCPM, color = Annotation)) +
    geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
    geom_smooth(method = "gam",  #loess for less than 1000 obs, gam for the more than 1000     
                aes(x = dptval, y = logCPM),
                color = "black",        # Smooth curve color
                size = 1.2,            # Curve line thickness
                se = TRUE) +          # Disable confidence intervals
    labs(
      title = paste0(GENE),
      x = "Pseudotime",
      y = "Normalized Expression (LogCPM)"
    ) +
    theme(
		panel.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(color = "black"), 
		axis.ticks = element_line(color = "black"),
		plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14), 
		legend.position = "none"
		) +
	scale_x_continuous(breaks = c(min(ConvertTable$dptval), max(ConvertTable$dptval)),
					   labels = c("Early", "Late")
						) +
	scale_color_manual(values = AnnotationColors,name = "Annotation")
  )
 }
dev.off()

	# Hh score vs. FHF pseudotime
UCell_score <- read.csv("./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

meta <- na.omit(FHF_Slingshot_SelectedCells@meta.data[, c("Annotation", "TimeStampMerged", "SlingshotPseudotime")])
colnames(meta) <- c("Annotation", "TimeStampMerged", "dptval")
meta_sort <- meta[order(meta$dptval),]

UCell_score_subset <- UCell_score[UCell_score$X %in% rownames(meta_sort), ]
UCell_score_subset <- UCell_score_subset[match(rownames(meta_sort), UCell_score_subset$X), ]
meta_sort$UCell <- UCell_score_subset$UCell

pdf(file = "./3-Merge3/Figures/Figure2_LimitedHhScore_FHF.pdf", width = 8, height = 8)
print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - FHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late"))
)
 print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 2.5, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - FHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")) +
  scale_y_continuous(limits = c(0, 0.25))
  )
p <- ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - FHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")
  )
print(p)

plot_data <- ggplot_build(p)$data[[1]]  # This extracts the first layer (the smooth curve)
write.csv(as.data.frame(plot_data), file = "./3-Merge3/Figures/Hh-Score_plot_data_FHF.csv")

print(ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.5) +  # Confidence interval ribbon
  geom_line(color = "black", size = 1.2) +  # Smooth curve
  labs(
    title = "Limited Hh score - FHF trajectory",
    x = "Pseudotime",
    y = "Smooth average Hedgehog Score (AU)"
  ) +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_x_continuous(
    breaks = c(min(plot_data$x), max(plot_data$x)),
    labels = c("Early", "Late")
  ) +
  scale_y_continuous(limits = c(0, 0.06))
)
dev.off()


#########################################
### JCF Slingshot trajectory analysis ###
#########################################
JCF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/JCF1/JCF_Slingshot_SelectedCells.RDS")

JCF_Slingshot_SelectedCells@meta.data <- JCF_Slingshot_SelectedCells@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters %in% c("0","15","17") ~ "Primitive streak",
	seurat_clusters %in% c("3","10","16") ~ "Nascent mesoderm",
	seurat_clusters %in% c("1","14") ~ "Mixed mesoderm",
	seurat_clusters %in% "13" ~ "Paraxial/Pharyngeal mesoderm",
	seurat_clusters %in% c("5","8") ~ "ExE mesoderm",
	seurat_clusters %in% "6" ~ "aSHF/Pharyngeal mesoderm",
	seurat_clusters %in% "9" ~ "pSHF/Pharyngeal mesoderm",
	seurat_clusters %in% c("2","11") ~ "Mesenchyme",
	seurat_clusters %in% "4" ~ "JCF/Mesenchyme",
	seurat_clusters %in% "12" ~ "Differentiating CMs",
	seurat_clusters %in% "7" ~ "CMs",
	TRUE ~ seurat_clusters
  ))
JCF_Slingshot_SelectedCells@meta.data$Annotation <- factor(
  JCF_Slingshot_SelectedCells@meta.data$Annotation,
  levels = c("Primitive streak", "Nascent mesoderm","Mixed mesoderm","Paraxial/Pharyngeal mesoderm","ExE mesoderm","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm","Mesenchyme","JCF/Mesenchyme","Differentiating CMs","CMs")
)

# UMAP selected cells
JCF1_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
JCF1_dataSub <- subset(x = JCF1_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,8,11), invert=FALSE)

good <- data.frame(cells = rownames(JCF1_dataSub@meta.data), annot = JCF1_dataSub@meta.data$Annotation)
all <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), annot = "Other")
rest <- all[!(all$cells %in% good$cells),]
full <- rbind(good,rest)

temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0

temp@meta.data$CustomColors <- full$annot[match(rownames(temp@meta.data), full$cells)]
AnnotColors <- c(AnnotationColors, "Other" = "#E4E4E4")

pdf(file="./3-Merge3/Figures/Figure2_UMAP_selectedCells_ForJCFtrajectory.pdf", width=10, height=10)
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
dev.off()

# UMAP plot pseudotime
pdf(file = "./3-Merge3/Figures/Figure2_TrajectoryViz_JCF.pdf", width = 10, height = 10)
  good2 <- data.frame(cells = colnames(JCF_Slingshot_SelectedCells), Pseudotime = JCF_Slingshot_SelectedCells@meta.data$SlingshotPseudotime)
  all2 <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), Pseudotime = NA)
  rest2 <- all2[!(all2$cells %in% good2$cells),]
  full2 <- rbind(good2,rest2)
temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
temp@meta.data$Pseudo <- full2$Pseudotime[match(rownames(temp@meta.data), full2$cells)]
print(FeaturePlot(temp, features = "Pseudo", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
dev.off()

#pdf(file="./3-Merge3/Figures/Figure2_TrajectoryViz_JCF.pdf", width=10, height=10)
#print(DimPlot(JCF_Slingshot_SelectedCells, reduction = "umap", group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(JCF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 2), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(JCF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(JCF_Slingshot_SelectedCells, reduction = "harmony",dims = c(2, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(FeaturePlot(JCF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "umap", order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#print(FeaturePlot(JCF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "harmony",dims = c(1, 2), order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#dev.off()

# Ridgeplot timestamp contribution to pseudotime
preTable <- JCF_Slingshot_SelectedCells@meta.data[,c("TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(preTable) <- c("TimeStampMerged","Pseudotime","Annotation")
pdf(file = "./3-Merge3/Figures/Figure2_JCF_RidgePlot.pdf", width = 8, height = 8)
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 1.5, size = 0.3) +
    scale_fill_gradientn(
      colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
      name = "Pseudotime"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove grid lines
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),  # Add space for arrow and title
      axis.text.y = element_text(size = 20, margin = margin(l = 10)),  # Push timepoint labels to the right
      axis.text.x = element_text(size = 20),  # Adjust size of x-axis labels
      axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1)  # Adjust x-axis title
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)),  # First and last breaks
      labels = c("Early", "Late")  # Replace numerical x-axis labels
    ) +
    labs(x = "Pseudotime", y = "Cell density") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    ) 
 )
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges(aes(fill = Annotation), scale = 1.5, size = 0.3) +
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
# Calculate densities and interpolate to match original data length
density_data <- preTable %>%
  group_by(TimeStampMerged) %>%
  mutate(
    Density = approxfun(density(Pseudotime, adjust = 1)$x, 
                        density(Pseudotime, adjust = 1)$y)(Pseudotime)
  )

# Plot with accurate density and coloring by Annotation
print(ggplot(density_data, aes(x = Pseudotime, y = TimeStampMerged, height = Density, fill = Annotation)) +
    geom_density_ridges(stat = "identity", scale = 1, size = 0.3) + # Use stat = "identity"
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
dev.off()

# Heatmap + GO terms top 500 DEGs along the pseudotime 
JCF_DEGlist <- read.csv("./3-Merge3/Slingshot/JCF1/tradeSeq/JCF_DEGlist.csv")
rownames(JCF_DEGlist) <- JCF_DEGlist$X
JCF_DEGlist <- JCF_DEGlist[order(-JCF_DEGlist$waldStat), ]
top500_JCF <- rownames(JCF_DEGlist[1:500,])

magict4_data <- readRDS(file = "./3-Merge3/DATANormcountsAfterMAGIC_t4.RDS")
meta <- JCF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
magict4_data <- magict4_data[, rownames(meta_sort)]
magict4_data <- magict4_data[top500_JCF, ]
tmagict4_data <- t(scale(t(magict4_data)))  # Row scaling (genes)

	#kmeans
set.seed(123)
wss <- sapply(1:15, 
              function(k){kmeans(tmagict4_data, k, nstart=25,iter.max = 15)$tot.withinss})
pdf(file = "./3-Merge3/Slingshot/JCF1/tradeSeq/JCF_HeatmapTop500_kmean_ElbowPlot.pdf", width = 8, height = 8)	
plot(1:15, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()
	
km_clusters <- data.frame(cluster = kmeans(tmagict4_data, centers = 6, iter.max = 15, nstart = 25)$cluster) # 25 random sets tested, 15 iteration for convergence
km_clusters$geneid <- rownames(km_clusters)
write.csv(as.data.frame(km_clusters), file = "./3-Merge3/Slingshot/JCF1/tradeSeq/Top500DEGs_km_clusters_JCF.csv")

row_order = hclust(dist(tmagict4_data))

	#GO terms
universe <- rownames(JCF_Slingshot_SelectedCells)
JCF_go <- compareCluster(geneid~cluster, data = km_clusters, fun = 'enrichGO', OrgDb  = "org.Mm.eg.db", ont='BP', qvalueCutoff = 0.05,  keyType = "SYMBOL", universe = universe) #readable = T,
write.csv(as.data.frame(JCF_go@compareClusterResult), file = "./3-Merge3/Slingshot/JCF1/tradeSeq/Top500DEGs_GOterms_JCF.csv")
		# extract top 5 terms per cluster
n <- 5
JCF_go_list <- JCF_go@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(qvalue, n = n) %>%  # Select top N based on qvalue
  summarise(GO_Terms = list(Description), .groups = "drop") %>%
  deframe()

	#Heatmap
		# Annotate go terms on the heatmap
		# Create annotation object using anno_textbox()
go_anno <- rowAnnotation(
  textbox = anno_textbox(
    km_clusters$cluster, 
    JCF_go_list, 
    #by = 'anno_block',
    word_wrap = T, 
    add_new_line = T,
   background_gp = gpar(fill = 'transparent')
  ),
  width = unit(4.5, "in")
)
combined_annotation <- HeatmapAnnotation(
  pseudotime = meta_sort$dptval,
  Annotation = meta_sort$Annotation,
  TimeStampMerged = meta_sort$TimeStampMerged,
  col = list(
    pseudotime = colorRamp2(
      seq(min(meta_sort$dptval), max(meta_sort$dptval), length.out = 256),  # Divide pseudotime into 256 gradient steps
      viridis(256, option = "plasma")  # Full plasma palette with 256 colors
    ),
    Annotation = AnnotationColors,
    TimeStampMerged = TimepointColors
  ),
  annotation_legend_param = list(
    pseudotime = list(title = "Pseudotime", at = c(min(meta_sort$dptval), max(meta_sort$dptval)),labels = c("Early", "Late")),
    Annotation = list(title = "Cell Type"),
    TimeStampMerged = list(title = "Timepoints")
  ),
  annotation_label = list(pseudotime = "Pseudotime", Annotation = "Cell type", TimeStampMerged = "Timepoints")
)
pdf(file = "./3-Merge3/Figures/Figure2_JCF_HeatmapTop500.pdf", width = 26, height = 8)	   
draw(Heatmap(
  tmagict4_data,
  cluster_columns = FALSE,  # No clustering for cells
  cluster_rows = TRUE,      # Cluster genes
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE, 
  show_column_names = FALSE,
  col = colorRamp2(c(-1.5,0,1.5), hcl_palette = "viridis"),
  row_split = km_clusters$cluster,
  row_gap = unit(1, "mm"),
  row_title = NULL,
  right_annotation = go_anno,
  heatmap_legend_param = list(
    title = "Expression levels",
    at = c(-2, 2),
    labels = c("Low", "High")),
    width = unit(16.6, "in")
),
annotation_legend_side = "right",
heatmap_legend_side = "right"
)
dev.off()

## Plot unique gene along Pseudotime ## Curve version2 =========
	# logCPM
rna_data <- JCF_Slingshot_SelectedCells[["RNA"]]@data
meta <- JCF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
rna_data <- as.matrix(rna_data[, rownames(meta_sort)])

GENE_LIST <- c("Mesp1","Nkx2-5","Hand1","Foxf1","Mab21l2")

pdf(file = "./3-Merge3/Figures/Figure2_JCF_GeneCurve.pdf", width = 8, height = 8)
for(GENE in GENE_LIST){
ConvertTable <- NULL
ConvertTable <- data.frame(logCPM = rna_data[rownames(rna_data) == GENE,], row.names = rownames(t(rna_data)))
ConvertTable$dptval <- meta_sort$dptval[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$Annotation <- meta_sort$Annotation[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$TimeStampMerged <- meta_sort$TimeStampMerged[match(rownames(ConvertTable), rownames(meta_sort))]
    print(ggplot(ConvertTable, aes(x = dptval, y = logCPM, color = Annotation)) +
    geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
    geom_smooth(method = "gam",  #loess for less than 1000 obs, gam for the more than 1000     
                aes(x = dptval, y = logCPM),
                color = "black",        # Smooth curve color
                size = 1.2,            # Curve line thickness
                se = TRUE) +          # Disable confidence intervals
    labs(
      title = paste0(GENE),
      x = "Pseudotime",
      y = "Normalized Expression (LogCPM)"
    ) +
    theme(
		panel.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(color = "black"), 
		axis.ticks = element_line(color = "black"),
		plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14), 
		legend.position = "none"
		) +
	scale_x_continuous(breaks = c(min(ConvertTable$dptval), max(ConvertTable$dptval)),
					   labels = c("Early", "Late")
						) +
	scale_color_manual(values = AnnotationColors,name = "Annotation")
  )
 }
dev.off()

	# Hh score vs. JCF pseudotime
UCell_score <- read.csv("./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

meta <- na.omit(JCF_Slingshot_SelectedCells@meta.data[, c("Annotation", "TimeStampMerged", "SlingshotPseudotime")])
colnames(meta) <- c("Annotation", "TimeStampMerged", "dptval")
meta_sort <- meta[order(meta$dptval),]

UCell_score_subset <- UCell_score[UCell_score$X %in% rownames(meta_sort), ]
UCell_score_subset <- UCell_score_subset[match(rownames(meta_sort), UCell_score_subset$X), ]
meta_sort$UCell <- UCell_score_subset$UCell

pdf(file = "./3-Merge3/Figures/Figure2_LimitedHhScore_JCF.pdf", width = 8, height = 8)
print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - JCF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late"))
)
 print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 2.5, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - JCF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")) +
  scale_y_continuous(limits = c(0, 0.25))
  )
p <- ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - JCF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")
  )
print(p)

plot_data <- ggplot_build(p)$data[[1]]  # This extracts the first layer (the smooth curve)
write.csv(as.data.frame(plot_data), file = "./3-Merge3/Figures/Hh-Score_plot_data_JCF.csv")

print(ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.5) +  # Confidence interval ribbon
  geom_line(color = "black", size = 1.2) +  # Smooth curve
  labs(
    title = "Limited Hh score - JCF trajectory",
    x = "Pseudotime",
    y = "Smooth average Hedgehog Score (AU)"
  ) +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_x_continuous(
    breaks = c(min(plot_data$x), max(plot_data$x)),
    labels = c("Early", "Late")
  ) +
  scale_y_continuous(limits = c(0, 0.06))
)
dev.off()

#########################################
### aSHF Slingshot trajectory analysis ###
#########################################
aSHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/SHF2/aSHF_Slingshot_SelectedCells.RDS")

aSHF_Slingshot_SelectedCells@meta.data <- aSHF_Slingshot_SelectedCells@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters %in% c("0","15","17") ~ "Primitive streak",
	seurat_clusters %in% c("3","10","16") ~ "Nascent mesoderm",
	seurat_clusters %in% c("1","14") ~ "Mixed mesoderm",
	seurat_clusters %in% "13" ~ "Paraxial/Pharyngeal mesoderm",
	seurat_clusters %in% c("5","8") ~ "ExE mesoderm",
	seurat_clusters %in% "6" ~ "aSHF/Pharyngeal mesoderm",
	seurat_clusters %in% "9" ~ "pSHF/Pharyngeal mesoderm",
	seurat_clusters %in% c("2","11") ~ "Mesenchyme",
	seurat_clusters %in% "4" ~ "JCF/Mesenchyme",
	seurat_clusters %in% "12" ~ "Differentiating CMs",
	seurat_clusters %in% "7" ~ "CMs",
	TRUE ~ seurat_clusters
  ))
aSHF_Slingshot_SelectedCells@meta.data$Annotation <- factor(
  aSHF_Slingshot_SelectedCells@meta.data$Annotation,
  levels = c("Primitive streak", "Nascent mesoderm","Mixed mesoderm","Paraxial/Pharyngeal mesoderm","ExE mesoderm","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm","Mesenchyme","JCF/Mesenchyme","Differentiating CMs","CMs")
)

# UMAP selected cells
SHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
SHF2_dataSub <- subset(x = SHF2_dataPreSub, subset = seurat_clusters %in% c(1,6,9), invert=FALSE)

good <- data.frame(cells = rownames(SHF2_dataSub@meta.data), annot = SHF2_dataSub@meta.data$Annotation)
all <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), annot = "Other")
rest <- all[!(all$cells %in% good$cells),]
full <- rbind(good,rest)

temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0

temp@meta.data$CustomColors <- full$annot[match(rownames(temp@meta.data), full$cells)]
AnnotColors <- c(AnnotationColors, "Other" = "#E4E4E4")

pdf(file="./3-Merge3/Figures/Figure2_UMAP_selectedCells_ForaSHFtrajectory.pdf", width=10, height=10)
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
dev.off()

# UMAP plot pseudotime
pdf(file = "./3-Merge3/Figures/Figure2_TrajectoryViz_aSHF.pdf", width = 10, height = 10)
  good2 <- data.frame(cells = colnames(aSHF_Slingshot_SelectedCells), Pseudotime = aSHF_Slingshot_SelectedCells@meta.data$SlingshotPseudotime)
  all2 <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), Pseudotime = NA)
  rest2 <- all2[!(all2$cells %in% good2$cells),]
  full2 <- rbind(good2,rest2)
temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
temp@meta.data$Pseudo <- full2$Pseudotime[match(rownames(temp@meta.data), full2$cells)]
print(FeaturePlot(temp, features = "Pseudo", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
dev.off()

#pdf(file="./3-Merge3/Figures/Figure2_TrajectoryViz_aSHF.pdf", width=10, height=10)
#print(DimPlot(aSHF_Slingshot_SelectedCells, reduction = "umap", group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(aSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 2), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(aSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(aSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(2, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(FeaturePlot(aSHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "umap", order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#print(FeaturePlot(aSHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "harmony",dims = c(1, 2), order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#dev.off()

# Ridgeplot timestamp contribution to pseudotime
preTable <- aSHF_Slingshot_SelectedCells@meta.data[,c("TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(preTable) <- c("TimeStampMerged","Pseudotime","Annotation")
pdf(file = "./3-Merge3/Figures/Figure2_aSHF_RidgePlot.pdf", width = 8, height = 8)
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 1.5, size = 0.3) +
    scale_fill_gradientn(
      colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
      name = "Pseudotime"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove grid lines
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),  # Add space for arrow and title
      axis.text.y = element_text(size = 20, margin = margin(l = 10)),  # Push timepoint labels to the right
      axis.text.x = element_text(size = 20),  # Adjust size of x-axis labels
      axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1)  # Adjust x-axis title
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)),  # First and last breaks
      labels = c("Early", "Late")  # Replace numerical x-axis labels
    ) +
    labs(x = "Pseudotime", y = "Cell density") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    ) 
 )
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges(aes(fill = Annotation), scale = 1.5, size = 0.3) +
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
# Calculate densities and interpolate to match original data length
density_data <- preTable %>%
  group_by(TimeStampMerged) %>%
  mutate(
    Density = approxfun(density(Pseudotime, adjust = 1)$x, 
                        density(Pseudotime, adjust = 1)$y)(Pseudotime)
  )

# Plot with accurate density and coloring by Annotation
print(ggplot(density_data, aes(x = Pseudotime, y = TimeStampMerged, height = Density, fill = Annotation)) +
    geom_density_ridges(stat = "identity", scale = 1, size = 0.3) + # Use stat = "identity"
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
dev.off()

# Heatmap + GO terms top 500 DEGs along the pseudotime 
aSHF_DEGlist <- read.csv("./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/aSHF_DEGlist.csv")
rownames(aSHF_DEGlist) <- aSHF_DEGlist$X
aSHF_DEGlist <- aSHF_DEGlist[order(-aSHF_DEGlist$waldStat), ]
top500_aSHF <- rownames(aSHF_DEGlist[1:500,])

magict4_data <- readRDS(file = "./3-Merge3/DATANormcountsAfterMAGIC_t4.RDS")
meta <- aSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
magict4_data <- magict4_data[, rownames(meta_sort)]
magict4_data <- magict4_data[top500_aSHF, ]
tmagict4_data <- t(scale(t(magict4_data)))  # Row scaling (genes)

	#kmeans
set.seed(123)
wss <- sapply(1:15, 
              function(k){kmeans(tmagict4_data, k, nstart=25,iter.max = 15)$tot.withinss})
pdf(file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/aSHF_HeatmapTop500_kmean_ElbowPlot.pdf", width = 8, height = 8)	
plot(1:15, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()
	
km_clusters <- data.frame(cluster = kmeans(tmagict4_data, centers = 6, iter.max = 15, nstart = 25)$cluster) # 25 random sets tested, 15 iteration for convergence
km_clusters$geneid <- rownames(km_clusters)
write.csv(as.data.frame(km_clusters), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/Top500DEGs_km_clusters_aSHF.csv")

row_order = hclust(dist(tmagict4_data))

	#GO terms
universe <- rownames(aSHF_Slingshot_SelectedCells)
aSHF_go <- compareCluster(geneid~cluster, data = km_clusters, fun = 'enrichGO', OrgDb  = "org.Mm.eg.db", ont='BP', qvalueCutoff = 0.05,  keyType = "SYMBOL", universe = universe) #readable = T,
write.csv(as.data.frame(aSHF_go@compareClusterResult), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/aSHF/Top500DEGs_GOterms_aSHF.csv")
		# extract top 5 terms per cluster
n <- 5
aSHF_go_list <- aSHF_go@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(qvalue, n = n) %>%  # Select top N based on qvalue
  summarise(GO_Terms = list(Description), .groups = "drop") %>%
  deframe()

	#Heatmap
		# Annotate go terms on the heatmap
		# Create annotation object using anno_textbox()
go_anno <- rowAnnotation(
  textbox = anno_textbox(
    km_clusters$cluster, 
    aSHF_go_list, 
    #by = 'anno_block',
    word_wrap = T, 
    add_new_line = T,
   background_gp = gpar(fill = 'transparent')
  ),
  width = unit(4.5, "in")
)
combined_annotation <- HeatmapAnnotation(
  pseudotime = meta_sort$dptval,
  Annotation = meta_sort$Annotation,
  TimeStampMerged = meta_sort$TimeStampMerged,
  col = list(
    pseudotime = colorRamp2(
      seq(min(meta_sort$dptval), max(meta_sort$dptval), length.out = 256),  # Divide pseudotime into 256 gradient steps
      viridis(256, option = "plasma")  # Full plasma palette with 256 colors
    ),
    Annotation = AnnotationColors,
    TimeStampMerged = TimepointColors
  ),
  annotation_legend_param = list(
    pseudotime = list(title = "Pseudotime", at = c(min(meta_sort$dptval), max(meta_sort$dptval)),labels = c("Early", "Late")),
    Annotation = list(title = "Cell Type"),
    TimeStampMerged = list(title = "Timepoints")
  ),
  annotation_label = list(pseudotime = "Pseudotime", Annotation = "Cell type", TimeStampMerged = "Timepoints")
)
pdf(file = "./3-Merge3/Figures/Figure2_aSHF_HeatmapTop500.pdf", width = 26, height = 8)	   
draw(Heatmap(
  tmagict4_data,
  cluster_columns = FALSE,  # No clustering for cells
  cluster_rows = TRUE,      # Cluster genes
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE, 
  show_column_names = FALSE,
  col = colorRamp2(c(-1.5,0,1.5), hcl_palette = "viridis"),
  row_split = km_clusters$cluster,
  row_gap = unit(1, "mm"),
  row_title = NULL,
  right_annotation = go_anno,
  heatmap_legend_param = list(
    title = "Expression levels",
    at = c(-2, 2),
    labels = c("Low", "High")),
    width = unit(16.6, "in")
),
annotation_legend_side = "right",
heatmap_legend_side = "right"
)
dev.off()

## Plot unique gene along Pseudotime ## Curve version2 =========
	# logCPM
rna_data <- aSHF_Slingshot_SelectedCells[["RNA"]]@data
meta <- aSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
rna_data <- as.matrix(rna_data[, rownames(meta_sort)])

GENE_LIST <- c("Isl1","Tbx1","Fgf8","Osr1","Aldh1a2")

pdf(file = "./3-Merge3/Figures/Figure2_aSHF_GeneCurve.pdf", width = 8, height = 8)
for(GENE in GENE_LIST){
ConvertTable <- NULL
ConvertTable <- data.frame(logCPM = rna_data[rownames(rna_data) == GENE,], row.names = rownames(t(rna_data)))
ConvertTable$dptval <- meta_sort$dptval[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$Annotation <- meta_sort$Annotation[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$TimeStampMerged <- meta_sort$TimeStampMerged[match(rownames(ConvertTable), rownames(meta_sort))]
    print(ggplot(ConvertTable, aes(x = dptval, y = logCPM, color = Annotation)) +
    geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
    geom_smooth(method = "gam",  #loess for less than 1000 obs, gam for the more than 1000     
                aes(x = dptval, y = logCPM),
                color = "black",        # Smooth curve color
                size = 1.2,            # Curve line thickness
                se = TRUE) +          # Disable confidence intervals
    labs(
      title = paste0(GENE),
      x = "Pseudotime",
      y = "Normalized Expression (LogCPM)"
    ) +
    theme(
		panel.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(color = "black"), 
		axis.ticks = element_line(color = "black"),
		plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14), 
		legend.position = "none"
		) +
	scale_x_continuous(breaks = c(min(ConvertTable$dptval), max(ConvertTable$dptval)),
					   labels = c("Early", "Late")
						) +
	scale_color_manual(values = AnnotationColors,name = "Annotation")
  )
 }
dev.off()

	# Hh score vs. aSHF pseudotime
UCell_score <- read.csv("./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

meta <- na.omit(aSHF_Slingshot_SelectedCells@meta.data[, c("Annotation", "TimeStampMerged", "SlingshotPseudotime")])
colnames(meta) <- c("Annotation", "TimeStampMerged", "dptval")
meta_sort <- meta[order(meta$dptval),]

UCell_score_subset <- UCell_score[UCell_score$X %in% rownames(meta_sort), ]
UCell_score_subset <- UCell_score_subset[match(rownames(meta_sort), UCell_score_subset$X), ]
meta_sort$UCell <- UCell_score_subset$UCell

pdf(file = "./3-Merge3/Figures/Figure2_LimitedHhScore_aSHF.pdf", width = 8, height = 8)
print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - aSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late"))
)
 print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 2.5, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - aSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")) +
  scale_y_continuous(limits = c(0, 0.25))
  )
p <- ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - aSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")
  )
print(p)

plot_data <- ggplot_build(p)$data[[1]]  # This extracts the first layer (the smooth curve)
write.csv(as.data.frame(plot_data), file = "./3-Merge3/Figures/Hh-Score_plot_data_aSHF.csv")

print(ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.5) +  # Confidence interval ribbon
  geom_line(color = "black", size = 1.2) +  # Smooth curve
  labs(
    title = "Limited Hh score - aSHF trajectory",
    x = "Pseudotime",
    y = "Smooth average Hedgehog Score (AU)"
  ) +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_x_continuous(
    breaks = c(min(plot_data$x), max(plot_data$x)),
    labels = c("Early", "Late")
  ) +
  scale_y_continuous(limits = c(0, 0.06))
)
dev.off()


#########################################
### pSHF Slingshot trajectory analysis ###
#########################################
pSHF_Slingshot_SelectedCells <- readRDS(file = "./3-Merge3/Slingshot/SHF2/pSHF_Slingshot_SelectedCells.RDS")

pSHF_Slingshot_SelectedCells@meta.data <- pSHF_Slingshot_SelectedCells@meta.data %>%
  mutate(Annotation = case_when(
    seurat_clusters %in% c("0","15","17") ~ "Primitive streak",
	seurat_clusters %in% c("3","10","16") ~ "Nascent mesoderm",
	seurat_clusters %in% c("1","14") ~ "Mixed mesoderm",
	seurat_clusters %in% "13" ~ "Paraxial/Pharyngeal mesoderm",
	seurat_clusters %in% c("5","8") ~ "ExE mesoderm",
	seurat_clusters %in% "6" ~ "aSHF/Pharyngeal mesoderm",
	seurat_clusters %in% "9" ~ "pSHF/Pharyngeal mesoderm",
	seurat_clusters %in% c("2","11") ~ "Mesenchyme",
	seurat_clusters %in% "4" ~ "JCF/Mesenchyme",
	seurat_clusters %in% "12" ~ "Differentiating CMs",
	seurat_clusters %in% "7" ~ "CMs",
	TRUE ~ seurat_clusters
  ))
pSHF_Slingshot_SelectedCells@meta.data$Annotation <- factor(
  pSHF_Slingshot_SelectedCells@meta.data$Annotation,
  levels = c("Primitive streak", "Nascent mesoderm","Mixed mesoderm","Paraxial/Pharyngeal mesoderm","ExE mesoderm","aSHF/Pharyngeal mesoderm","pSHF/Pharyngeal mesoderm","Mesenchyme","JCF/Mesenchyme","Differentiating CMs","CMs")
)

# UMAP selected cells
SHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
SHF2_dataSub <- subset(x = SHF2_dataPreSub, subset = seurat_clusters %in% c(1,6,9), invert=FALSE)

good <- data.frame(cells = rownames(SHF2_dataSub@meta.data), annot = SHF2_dataSub@meta.data$Annotation)
all <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), annot = "Other")
rest <- all[!(all$cells %in% good$cells),]
full <- rbind(good,rest)

temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0

temp@meta.data$CustomColors <- full$annot[match(rownames(temp@meta.data), full$cells)]
AnnotColors <- c(AnnotationColors, "Other" = "#E4E4E4")

pdf(file="./3-Merge3/Figures/Figure2_UMAP_selectedCells_ForpSHFtrajectory.pdf", width=10, height=10)
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors))
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = T, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
print(DimPlot(temp, reduction = "umap", group.by = "CustomColors", label = F, order = setdiff(names(AnnotColors), "Other"), pt.size = 0.7,raster = FALSE, cols = AnnotColors) + NoLegend())
dev.off()

# UMAP plot pseudotime
pdf(file = "./3-Merge3/Figures/Figure2_TrajectoryViz_pSHF.pdf", width = 10, height = 10)
  good2 <- data.frame(cells = colnames(pSHF_Slingshot_SelectedCells), Pseudotime = pSHF_Slingshot_SelectedCells@meta.data$SlingshotPseudotime)
  all2 <- data.frame(cells = rownames(Merge3_Krup_Tyser_Pijuan_harmony2_res1.0@meta.data), Pseudotime = NA)
  rest2 <- all2[!(all2$cells %in% good2$cells),]
  full2 <- rbind(good2,rest2)
temp <- Merge3_Krup_Tyser_Pijuan_harmony2_res1.0
temp@meta.data$Pseudo <- full2$Pseudotime[match(rownames(temp@meta.data), full2$cells)]
print(FeaturePlot(temp, features = "Pseudo", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
dev.off()

#pdf(file="./3-Merge3/Figures/Figure2_TrajectoryViz_pSHF.pdf", width=10, height=10)
#print(DimPlot(pSHF_Slingshot_SelectedCells, reduction = "umap", group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(pSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 2), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(pSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(1, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(DimPlot(pSHF_Slingshot_SelectedCells, reduction = "harmony",dims = c(2, 3), group.by = "Annotation", pt.size = 0.7,raster = FALSE, cols = AnnotationColors))
#print(FeaturePlot(pSHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "umap", order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#print(FeaturePlot(pSHF_Slingshot_SelectedCells, features = "SlingshotPseudotime", reduction = "harmony",dims = c(1, 2), order = TRUE, cols = CustomPalette(low = "#0D0887FF", high = "#F0F921FF", mid = "#CC4678FF", k = 50)))
#dev.off()

# Ridgeplot timestamp contribution to pseudotime
preTable <- pSHF_Slingshot_SelectedCells@meta.data[,c("TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(preTable) <- c("TimeStampMerged","Pseudotime","Annotation")
pdf(file = "./3-Merge3/Figures/Figure2_pSHF_RidgePlot.pdf", width = 8, height = 8)
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges_gradient(aes(fill = ..x..), scale = 1.5, size = 0.3) +
    scale_fill_gradientn(
      colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF"),
      name = "Pseudotime"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove grid lines
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),  # Add space for arrow and title
      axis.text.y = element_text(size = 20, margin = margin(l = 10)),  # Push timepoint labels to the right
      axis.text.x = element_text(size = 20),  # Adjust size of x-axis labels
      axis.title.x = element_text(size = 20, hjust = 0.5, vjust = 1)  # Adjust x-axis title
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)),  # First and last breaks
      labels = c("Early", "Late")  # Replace numerical x-axis labels
    ) +
    labs(x = "Pseudotime", y = "Cell density") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    ) 
 )
print(ggplot(preTable, aes(x = Pseudotime, y = TimeStampMerged)) +
    geom_density_ridges(aes(fill = Annotation), scale = 1.5, size = 0.3) +
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
# Calculate densities and interpolate to match original data length
density_data <- preTable %>%
  group_by(TimeStampMerged) %>%
  mutate(
    Density = approxfun(density(Pseudotime, adjust = 1)$x, 
                        density(Pseudotime, adjust = 1)$y)(Pseudotime)
  )

# Plot with accurate density and coloring by Annotation
print(ggplot(density_data, aes(x = Pseudotime, y = TimeStampMerged, height = Density, fill = Annotation)) +
    geom_density_ridges(stat = "identity", scale = 1, size = 0.3) + # Use stat = "identity"
    scale_fill_manual(
      values = AnnotationColors,
      name = "Annotation"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, hjust = 0.5, vjust = 1)
    ) +
    scale_x_continuous(
      breaks = c(min(preTable$Pseudotime), max(preTable$Pseudotime)), 
      labels = c("Early", "Late")
    ) +
    labs(x = "Pseudotime") +
    theme(
      axis.line.x = element_line(
        arrow = arrow(type = "closed", length = unit(0.25, "cm"), ends = "last"),
        colour = "black"
      )
    )
)
dev.off()

# Heatmap + GO terms top 500 DEGs along the pseudotime 
pSHF_DEGlist <- read.csv("./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/pSHF_DEGlist.csv")
rownames(pSHF_DEGlist) <- pSHF_DEGlist$X
pSHF_DEGlist <- pSHF_DEGlist[order(-pSHF_DEGlist$waldStat), ]
top500_pSHF <- rownames(pSHF_DEGlist[1:500,])

magict4_data <- readRDS(file = "./3-Merge3/DATANormcountsAfterMAGIC_t4.RDS")
meta <- pSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
magict4_data <- magict4_data[, rownames(meta_sort)]
magict4_data <- magict4_data[top500_pSHF, ]
tmagict4_data <- t(scale(t(magict4_data)))  # Row scaling (genes)

	#kmeans
set.seed(123)
wss <- sapply(1:15, 
              function(k){kmeans(tmagict4_data, k, nstart=25,iter.max = 15)$tot.withinss})
pdf(file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/pSHF_HeatmapTop500_kmean_ElbowPlot.pdf", width = 8, height = 8)	
plot(1:15, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()
	
km_clusters <- data.frame(cluster = kmeans(tmagict4_data, centers = 6, iter.max = 15, nstart = 25)$cluster) # 25 random sets tested, 15 iteration for convergence
km_clusters$geneid <- rownames(km_clusters)
write.csv(as.data.frame(km_clusters), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/Top500DEGs_km_clusters_pSHF.csv")

row_order = hclust(dist(tmagict4_data))

	#GO terms
universe <- rownames(pSHF_Slingshot_SelectedCells)
pSHF_go <- compareCluster(geneid~cluster, data = km_clusters, fun = 'enrichGO', OrgDb  = "org.Mm.eg.db", ont='BP', qvalueCutoff = 0.05,  keyType = "SYMBOL", universe = universe) #readable = T,
write.csv(as.data.frame(pSHF_go@compareClusterResult), file = "./3-Merge3/Slingshot/SHF2/tradeSeq/pSHF/Top500DEGs_GOterms_pSHF.csv")
		# extract top 5 terms per cluster
n <- 5
pSHF_go_list <- pSHF_go@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(qvalue, n = n) %>%  # Select top N based on qvalue
  summarise(GO_Terms = list(Description), .groups = "drop") %>%
  deframe()

	#Heatmap
		# Annotate go terms on the heatmap
		# Create annotation object using anno_textbox()
go_anno <- rowAnnotation(
  textbox = anno_textbox(
    km_clusters$cluster, 
    pSHF_go_list, 
    #by = 'anno_block',
    word_wrap = T, 
    add_new_line = T,
   background_gp = gpar(fill = 'transparent')
  ),
  width = unit(4.5, "in")
)
combined_annotation <- HeatmapAnnotation(
  pseudotime = meta_sort$dptval,
  Annotation = meta_sort$Annotation,
  TimeStampMerged = meta_sort$TimeStampMerged,
  col = list(
    pseudotime = colorRamp2(
      seq(min(meta_sort$dptval), max(meta_sort$dptval), length.out = 256),  # Divide pseudotime into 256 gradient steps
      viridis(256, option = "plasma")  # Full plasma palette with 256 colors
    ),
    Annotation = AnnotationColors,
    TimeStampMerged = TimepointColors
  ),
  annotation_legend_param = list(
    pseudotime = list(title = "Pseudotime", at = c(min(meta_sort$dptval), max(meta_sort$dptval)),labels = c("Early", "Late")),
    Annotation = list(title = "Cell Type"),
    TimeStampMerged = list(title = "Timepoints")
  ),
  annotation_label = list(pseudotime = "Pseudotime", Annotation = "Cell type", TimeStampMerged = "Timepoints")
)
pdf(file = "./3-Merge3/Figures/Figure2_pSHF_HeatmapTop500.pdf", width = 26, height = 8)	   
draw(Heatmap(
  tmagict4_data,
  cluster_columns = FALSE,  # No clustering for cells
  cluster_rows = TRUE,      # Cluster genes
  top_annotation = combined_annotation,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE, 
  show_column_names = FALSE,
  col = colorRamp2(c(-1.5,0,1.5), hcl_palette = "viridis"),
  row_split = km_clusters$cluster,
  row_gap = unit(1, "mm"),
  row_title = NULL,
  right_annotation = go_anno,
  heatmap_legend_param = list(
    title = "Expression levels",
    at = c(-2, 2),
    labels = c("Low", "High")),
    width = unit(16.6, "in")
),
annotation_legend_side = "right",
heatmap_legend_side = "right"
)
dev.off()

## Plot unique gene along Pseudotime ## Curve version2 =========
	# logCPM
rna_data <- pSHF_Slingshot_SelectedCells[["RNA"]]@data
meta <- pSHF_Slingshot_SelectedCells@meta.data[,c("seurat_clusters","TimeStampMerged","SlingshotPseudotime","Annotation")]
colnames(meta) <- c("seurat_clusters","TimeStampMerged","dptval","Annotation")
meta_sort <- meta[order(meta$dptval),]
rna_data <- as.matrix(rna_data[, rownames(meta_sort)])

GENE_LIST <- c("Isl1","Tbx1","Fgf8","Osr1","Aldh1a2")

pdf(file = "./3-Merge3/Figures/Figure2_pSHF_GeneCurve.pdf", width = 8, height = 8)
for(GENE in GENE_LIST){
ConvertTable <- NULL
ConvertTable <- data.frame(logCPM = rna_data[rownames(rna_data) == GENE,], row.names = rownames(t(rna_data)))
ConvertTable$dptval <- meta_sort$dptval[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$Annotation <- meta_sort$Annotation[match(rownames(ConvertTable), rownames(meta_sort))]
ConvertTable$TimeStampMerged <- meta_sort$TimeStampMerged[match(rownames(ConvertTable), rownames(meta_sort))]
    print(ggplot(ConvertTable, aes(x = dptval, y = logCPM, color = Annotation)) +
    geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
    geom_smooth(method = "gam",  #loess for less than 1000 obs, gam for the more than 1000     
                aes(x = dptval, y = logCPM),
                color = "black",        # Smooth curve color
                size = 1.2,            # Curve line thickness
                se = TRUE) +          # Disable confidence intervals
    labs(
      title = paste0(GENE),
      x = "Pseudotime",
      y = "Normalized Expression (LogCPM)"
    ) +
    theme(
		panel.background = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(color = "black"), 
		axis.ticks = element_line(color = "black"),
		plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14), 
		legend.position = "none"
		) +
	scale_x_continuous(breaks = c(min(ConvertTable$dptval), max(ConvertTable$dptval)),
					   labels = c("Early", "Late")
						) +
	scale_color_manual(values = AnnotationColors,name = "Annotation")
  )
 }
dev.off()

	# Hh score vs. pSHF pseudotime
UCell_score <- read.csv("./3-Merge3/UCell/UCell_score_Merge3_Krup_Tyser_Pijuan_harmony2_res1.0_AfterPrepSCTFindMarkers.csv")

meta <- na.omit(pSHF_Slingshot_SelectedCells@meta.data[, c("Annotation", "TimeStampMerged", "SlingshotPseudotime")])
colnames(meta) <- c("Annotation", "TimeStampMerged", "dptval")
meta_sort <- meta[order(meta$dptval),]

UCell_score_subset <- UCell_score[UCell_score$X %in% rownames(meta_sort), ]
UCell_score_subset <- UCell_score_subset[match(rownames(meta_sort), UCell_score_subset$X), ]
meta_sort$UCell <- UCell_score_subset$UCell

pdf(file = "./3-Merge3/Figures/Figure2_LimitedHhScore_pSHF.pdf", width = 8, height = 8)
print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - pSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late"))
)
 print(ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_point(size = 4, alpha = 0.5) +  # Scatter plot points
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 2.5, se = TRUE) +  # Smooth curve
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - pSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")) +
  scale_y_continuous(limits = c(0, 0.25))
  )
p <- ggplot(meta_sort, aes(x = dptval, y = UCell, color = Annotation)) +
  geom_smooth(method = "gam", aes(x = dptval, y = UCell), color = "black", size = 1.2, se = TRUE) +
  scale_color_manual(values = AnnotationColors) +
  labs(
    title = "Limited Hh score - pSHF trajectory",
    x = "Pseudotime",
    y = "Hedgehog Score (AU)"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  scale_x_continuous(
    breaks = c(min(meta_sort$dptval), max(meta_sort$dptval)),
    labels = c("Early", "Late")
  )
print(p)

plot_data <- ggplot_build(p)$data[[1]]  # This extracts the first layer (the smooth curve)
write.csv(as.data.frame(plot_data), file = "./3-Merge3/Figures/Hh-Score_plot_data_pSHF.csv")

print(ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey80", alpha = 0.5) +  # Confidence interval ribbon
  geom_line(color = "black", size = 1.2) +  # Smooth curve
  labs(
    title = "Limited Hh score - pSHF trajectory",
    x = "Pseudotime",
    y = "Smooth average Hedgehog Score (AU)"
  ) +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_x_continuous(
    breaks = c(min(plot_data$x), max(plot_data$x)),
    labels = c("Early", "Late")
  ) +
  scale_y_continuous(limits = c(0, 0.06))
)
dev.off()


	# Hh score vs. all pseudotime
Score_plot_data_FHF <- read.csv("./3-Merge3/Figures/Hh-Score_plot_data_FHF.csv", row.names=1)
Score_plot_data_FHF$Annotation <- "FHF"
Score_plot_data_FHF$X <- seq(1,80)
Score_plot_data_JCF <- read.csv("./3-Merge3/Figures/Hh-Score_plot_data_JCF.csv", row.names=1)
Score_plot_data_JCF$Annotation <- "JCF"
Score_plot_data_JCF$X <- seq(1,80)
Score_plot_data_aSHF <- read.csv("./3-Merge3/Figures/Hh-Score_plot_data_aSHF.csv", row.names=1)
Score_plot_data_aSHF$Annotation <- "aSHF"
Score_plot_data_aSHF$X <- seq(1,80)
Score_plot_data_pSHF <- read.csv("./3-Merge3/Figures/Hh-Score_plot_data_pSHF.csv", row.names=1)
Score_plot_data_pSHF$Annotation <- "pSHF"
Score_plot_data_pSHF$X <- seq(1,80)

Merge_Score_plot_data <- rbind(Score_plot_data_FHF,Score_plot_data_JCF,Score_plot_data_aSHF,Score_plot_data_pSHF)

Merge_Score_plot_data$Annotation <- factor(
  Merge_Score_plot_data$Annotation,
  levels = c("FHF","aSHF","pSHF","JCF")
)

pdf(file = "./3-Merge3/Figures/Figure2_LimitedHhScore_All-Trajectories.pdf", width = 8, height = 8)
print(ggplot(Merge_Score_plot_data, aes(x = X, y = y, color = Annotation, group = Annotation)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, group = Annotation), fill = "grey80", alpha = 0.5) +  # Confidence interval ribbon
  geom_line(aes(group = Annotation, color = Annotation),size = 1.2) +  # Smooth curve
  scale_color_manual(name = "Trajectories", values = c("FHF"="#E41A1C","aSHF"="#40928a","pSHF"="#5e4fa2","JCF"="#fb9ac2")) +
  labs(
    title = "",
    x = "Pseudotime",
    y = "Smooth average Hedgehog Score (AU)"
  ) +
    theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_x_continuous(
    breaks = c(min(Merge_Score_plot_data$x), max(Merge_Score_plot_data$x)),
    labels = c("Early", "Late")
  ) +
  scale_y_continuous(limits = c(0, 0.06))
)
dev.off()











