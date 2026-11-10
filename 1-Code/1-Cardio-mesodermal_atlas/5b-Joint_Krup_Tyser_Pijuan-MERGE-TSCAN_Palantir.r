##################################################################################################################################################################################################
################################                                                                                                                                  ################################
################################ Krup-Pijuan-Tyser integrated Atlas - Merged dataset - TSCAN and Palantir trajectory analyses and comparison to Slingshot results ################################
################################                                                                                                                                  ################################
##################################################################################################################################################################################################
# R/4.3.1
library(Seurat)
library(TSCAN)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(dplyr)
library(TrajectoryUtils)

Merge3_Krup_Tyser_Pijuan_harmony2_res1.0 <- readRDS(file = "./3-Merge3/Merge3_Krup_Tyser_Pijuan_harmony2_res1.0.RDS")

#===============================================================================================================================================#
#####################################################========      TSCAN trajectory analysis      ========#######################################
#===============================================================================================================================================#

	#######
	# FHF #
	#######
FHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"), invert=FALSE)
FHF2_dataSub <- subset(x = FHF2_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,7,8,11,12), invert=FALSE)

# ─────────────────────────────────────────────
# 1. Extract Harmony embeddings from Seurat
# ─────────────────────────────────────────────
harmony_embeddings <- Embeddings(FHF2_dataSub, reduction = "harmony")[, 1:30]

# ─────────────────────────────────────────────
# 2. Build a SCE and inject Harmony + clusters
# ─────────────────────────────────────────────
sce <- as.SingleCellExperiment(FHF2_dataSub)

# Store the 30 Harmony PCs as a reducedDim named "HARMONY"
reducedDim(sce, "HARMONY") <- harmony_embeddings

# Pull Seurat cluster labels (adjust slot name if needed)
sce$seurat_clusters <- FHF2_dataSub$seurat_clusters   # or $annotated_clusters
sce$seurat_clusters <- droplevels(sce$seurat_clusters)

# ─────────────────────────────────────────────
# 3. Build the MST with TSCAN using Harmony space
# ─────────────────────────────────────────────
# Compute per-cluster centroids in Harmony space TSCAN does this internally when pass `use.dimred`
set.seed(42)
mst <- createClusterMST(
  sce,
  clusters   = sce$seurat_clusters,   # annotation, not TSCAN's own clustering
  use.dimred = "HARMONY"              # Harmony PCs
)

# Inspect the MST
print(mst)
pdf(file = "./TSCAN/1-FHF_mst.pdf", width = 5, height = 5)
plot(mst, main = "TSCAN MST (Harmony space, Seurat clusters)")
dev.off()

# ─────────────────────────────────────────────
# 4. Define paths through the MST
# ─────────────────────────────────────────────
root <- "1"
paths <- TrajectoryUtils::defineMSTPaths(
  g     = mst,
  roots = root
)

# Build cells × paths pseudotime matrix
harmony_mat <- reducedDim(sce, "HARMONY")
clusters    <- sce$seurat_clusters

centroids <- t(sapply(levels(clusters), function(cl) {
  colMeans(harmony_mat[clusters == cl, , drop = FALSE])
}))
rownames(centroids) <- levels(clusters)

pseudo_mat <- matrix(
  NA_real_,
  nrow     = ncol(sce),
  ncol     = length(paths),
  dimnames = list(colnames(sce), paste0("path", seq_along(paths)))
)

for (i in seq_along(paths)) {
  path_clusters  <- paths[[i]]
  path_centroids <- centroids[path_clusters, , drop = FALSE]
  cum_dist <- c(0, cumsum(sapply(2:nrow(path_centroids), function(j) {
    sqrt(sum((path_centroids[j, ] - path_centroids[j-1, ])^2))
  })))
  names(cum_dist) <- path_clusters
  
  for (j in seq_along(path_clusters)) {
    cl       <- path_clusters[j]
    cell_idx <- which(clusters == cl)
    if (length(cell_idx) == 0) next
    
    if (j < length(path_clusters)) {
      seg_start <- path_centroids[j, ]
      seg_end   <- path_centroids[j+1, ]
      seg_unit  <- (seg_end - seg_start) / sqrt(sum((seg_end - seg_start)^2))
      projections <- as.numeric((harmony_mat[cell_idx, , drop = FALSE] - 
                       rep(seg_start, each = length(cell_idx))) %*% seg_unit)
      projections <- pmax(0, pmin(projections, sqrt(sum((seg_end - seg_start)^2))))
      pseudo_mat[cell_idx, i] <- cum_dist[j] + projections
    } else {
      pseudo_mat[cell_idx, i] <- cum_dist[j]
    }
  }
}

# Store in SCE and Seurat object
for (i in seq_len(ncol(pseudo_mat))) {
  col_name               <- paste0("TSCAN_path", i)
  sce[[col_name]]        <- pseudo_mat[, i]
  FHF2_dataSub[[col_name]] <- pseudo_mat[, i]
}

saveRDS(FHF2_dataSub, file = "./TSCAN/1-FHF2_dataSub_TSCAN.RDS")

#R/4.0.3
#FHF2_dataSub_TSCAN <- readRDS(file = "./TSCAN/1-FHF2_dataSub_TSCAN.RDS")
#pdf(file="./TSCAN/1-FHF2_dataSub_TSCAN_trajectories.pdf", width=10, height=10)
#print(FeaturePlot(FHF2_dataSub_TSCAN, features = "TSCAN_path1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#print(FeaturePlot(FHF2_dataSub_TSCAN, features = "TSCAN_path2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#dev.off()


	#######
	# JCF #
	#######
JCF1_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
JCF1_dataSub <- subset(x = JCF1_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,8,11), invert=FALSE)

# ─────────────────────────────────────────────
# 1. Extract Harmony embeddings from Seurat
# ─────────────────────────────────────────────
harmony_embeddings <- Embeddings(JCF1_dataSub, reduction = "harmony")[, 1:30]

# ─────────────────────────────────────────────
# 2. Build a SCE and inject Harmony + clusters
# ─────────────────────────────────────────────
sce <- as.SingleCellExperiment(JCF1_dataSub)

# Store the 30 Harmony PCs as a reducedDim named "HARMONY"
reducedDim(sce, "HARMONY") <- harmony_embeddings

# Pull Seurat cluster labels (adjust slot name if needed)
sce$seurat_clusters <- JCF1_dataSub$seurat_clusters   # or $annotated_clusters
sce$seurat_clusters <- droplevels(sce$seurat_clusters)

# ─────────────────────────────────────────────
# 3. Build the MST with TSCAN using Harmony space
# ─────────────────────────────────────────────
# Compute per-cluster centroids in Harmony space TSCAN does this internally when pass `use.dimred`
set.seed(42)
mst <- createClusterMST(
  sce,
  clusters   = sce$seurat_clusters,   # annotation, not TSCAN's own clustering
  use.dimred = "HARMONY"              # requirement 1: Harmony PCs
)

# Inspect the MST
print(mst)
pdf(file = "./TSCAN/2-JCF_mst.pdf", width = 5, height = 5)
plot(mst, main = "TSCAN MST (Harmony space, Seurat clusters)")
dev.off()

# ─────────────────────────────────────────────
# 4. Define paths through the MST
# ─────────────────────────────────────────────
root <- "1"
paths <- TrajectoryUtils::defineMSTPaths(
  g     = mst,
  roots = root
)

# Build cells × paths pseudotime matrix
harmony_mat <- reducedDim(sce, "HARMONY")
clusters    <- sce$seurat_clusters

centroids <- t(sapply(levels(clusters), function(cl) {
  colMeans(harmony_mat[clusters == cl, , drop = FALSE])
}))
rownames(centroids) <- levels(clusters)

pseudo_mat <- matrix(
  NA_real_,
  nrow     = ncol(sce),
  ncol     = length(paths),
  dimnames = list(colnames(sce), paste0("path", seq_along(paths)))
)

for (i in seq_along(paths)) {
  path_clusters  <- paths[[i]]
  path_centroids <- centroids[path_clusters, , drop = FALSE]
  cum_dist <- c(0, cumsum(sapply(2:nrow(path_centroids), function(j) {
    sqrt(sum((path_centroids[j, ] - path_centroids[j-1, ])^2))
  })))
  names(cum_dist) <- path_clusters
  
  for (j in seq_along(path_clusters)) {
    cl       <- path_clusters[j]
    cell_idx <- which(clusters == cl)
    if (length(cell_idx) == 0) next
    
    if (j < length(path_clusters)) {
      seg_start <- path_centroids[j, ]
      seg_end   <- path_centroids[j+1, ]
      seg_unit  <- (seg_end - seg_start) / sqrt(sum((seg_end - seg_start)^2))
      projections <- as.numeric((harmony_mat[cell_idx, , drop = FALSE] - 
                       rep(seg_start, each = length(cell_idx))) %*% seg_unit)
      projections <- pmax(0, pmin(projections, sqrt(sum((seg_end - seg_start)^2))))
      pseudo_mat[cell_idx, i] <- cum_dist[j] + projections
    } else {
      pseudo_mat[cell_idx, i] <- cum_dist[j]
    }
  }
}

# Store in SCE and Seurat object
for (i in seq_len(ncol(pseudo_mat))) {
  col_name               <- paste0("TSCAN_path", i)
  sce[[col_name]]        <- pseudo_mat[, i]
  JCF1_dataSub[[col_name]] <- pseudo_mat[, i]
}

saveRDS(JCF1_dataSub, file = "./TSCAN/2-JCF1_dataSub_TSCAN.RDS")

#R/4.0.3
#JCF1_dataSub_TSCAN <- readRDS(file = "/project/imoskowitz/rcanac/12-SandwichModelPaper/3-Reviewing_Circulation/TSCAN/2-JCF1_dataSub_TSCAN.RDS")
#pdf(file="./TSCAN/2-JCF1_dataSub_TSCAN_trajectories.pdf", width=10, height=10)
#print(FeaturePlot(JCF1_dataSub_TSCAN, features = "TSCAN_path1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#print(FeaturePlot(JCF1_dataSub_TSCAN, features = "TSCAN_path2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#dev.off()


	#######
	# SHF #
	#######
SHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
SHF2_dataSub <- subset(x = SHF2_dataPreSub, subset = seurat_clusters %in% c(1,6,9), invert=FALSE)

# ─────────────────────────────────────────────
# 1. Extract Harmony embeddings from Seurat
# ─────────────────────────────────────────────
harmony_embeddings <- Embeddings(SHF2_dataSub, reduction = "harmony")[, 1:30]

# ─────────────────────────────────────────────
# 2. Build a SCE and inject Harmony + clusters
# ─────────────────────────────────────────────
sce <- as.SingleCellExperiment(SHF2_dataSub)

# Store the 30 Harmony PCs as a reducedDim named "HARMONY"
reducedDim(sce, "HARMONY") <- harmony_embeddings

# Pull Seurat cluster labels (adjust slot name if needed)
sce$seurat_clusters <- SHF2_dataSub$seurat_clusters   # or $annotated_clusters
sce$seurat_clusters <- droplevels(sce$seurat_clusters)

# ─────────────────────────────────────────────
# 3. Build the MST with TSCAN using Harmony space
# ─────────────────────────────────────────────
# Compute per-cluster centroids in Harmony space TSCAN does this internally when pass `use.dimred`
set.seed(42)
mst <- createClusterMST(
  sce,
  clusters   = sce$seurat_clusters,   # annotation, not TSCAN's own clustering
  use.dimred = "HARMONY"              # requirement 1: Harmony PCs
)

# Inspect the MST
print(mst)
pdf(file = "./TSCAN/3-SHF_mst.pdf", width = 5, height = 5)
plot(mst, main = "TSCAN MST (Harmony space, Seurat clusters)")
dev.off()

# ─────────────────────────────────────────────
# 4. Define paths through the MST
# ─────────────────────────────────────────────
paths <- list(
  path1 = c("1", "6"),
  path2 = c("1", "9")
)

# Build cells × paths pseudotime matrix
harmony_mat <- reducedDim(sce, "HARMONY")
clusters    <- sce$seurat_clusters

centroids <- t(sapply(levels(clusters), function(cl) {
  colMeans(harmony_mat[clusters == cl, , drop = FALSE])
}))
rownames(centroids) <- levels(clusters)

pseudo_mat <- matrix(
  NA_real_,
  nrow     = ncol(sce),
  ncol     = length(paths),
  dimnames = list(colnames(sce), paste0("path", seq_along(paths)))
)

for (i in seq_along(paths)) {
  path_clusters  <- paths[[i]]
  path_centroids <- centroids[path_clusters, , drop = FALSE]
  cum_dist <- c(0, cumsum(sapply(2:nrow(path_centroids), function(j) {
    sqrt(sum((path_centroids[j, ] - path_centroids[j-1, ])^2))
  })))
  names(cum_dist) <- path_clusters
  
  for (j in seq_along(path_clusters)) {
    cl       <- path_clusters[j]
    cell_idx <- which(clusters == cl)
    if (length(cell_idx) == 0) next
    
    if (j < length(path_clusters)) {
      seg_start <- path_centroids[j, ]
      seg_end   <- path_centroids[j+1, ]
      seg_unit  <- (seg_end - seg_start) / sqrt(sum((seg_end - seg_start)^2))
      projections <- as.numeric((harmony_mat[cell_idx, , drop = FALSE] - 
                       rep(seg_start, each = length(cell_idx))) %*% seg_unit)
      projections <- pmax(0, pmin(projections, sqrt(sum((seg_end - seg_start)^2))))
      pseudo_mat[cell_idx, i] <- cum_dist[j] + projections
    } else {
      pseudo_mat[cell_idx, i] <- cum_dist[j]
    }
  }
}

# Store in SCE and Seurat object
for (i in seq_len(ncol(pseudo_mat))) {
  col_name               <- paste0("TSCAN_path", i)
  sce[[col_name]]        <- pseudo_mat[, i]
  SHF2_dataSub[[col_name]] <- pseudo_mat[, i]
}

saveRDS(SHF2_dataSub, file = "./TSCAN/3-SHF2_dataSub_TSCAN.RDS")

#R/4.0.3
#SHF2_dataSub_TSCAN <- readRDS(file = "./TSCAN/3-SHF2_dataSub_TSCAN.RDS")
#pdf(file="./TSCAN/3-SHF2_dataSub_TSCAN_TSCAN_trajectories.pdf", width=10, height=10)
#print(FeaturePlot(SHF2_dataSub_TSCAN, features = "TSCAN_path1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#print(FeaturePlot(SHF2_dataSub_TSCAN, features = "TSCAN_path2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(color = "Pseudotime"))
#dev.off()

#==============================================================================================================================================#
###################################################========      Palantir trajectory analysis      ========#####################################
#==============================================================================================================================================#
library(reticulate)
use_condaenv("seuratextend", required = TRUE)
library(SeuratExtend)


	#######
	# FHF #
	#######
FHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E6.5","E6.75","E7.0","E7.25","E7.5"), invert=FALSE)
FHF2_dataSub <- subset(x = FHF2_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,7,8,11,12), invert=FALSE)

seurat_obj <- Palantir.RunDM(
  FHF2_dataSub,
  reduction = "harmony",   
  n_components = 30
)
seurat_obj@reductions

# Visualise the multiscale space
pdf(file="./Palantir/1-FHF2_dataSub_Palantir_ms-dm.pdf", width=10, height=10)
DimPlot2(seurat_obj, reduction = "ms", group.by = "seurat_clusters", label = TRUE)
DimPlot2(seurat_obj, reduction = "dm", group.by = "seurat_clusters", label = TRUE)
dev.off()

# Pick tip cell as root
ms_df <- as.data.frame(Embeddings(seurat_obj, "ms")[, 1:2])
colnames(ms_df) <- c("MS_1", "MS_2")
ms_df$cluster   <- seurat_obj$seurat_clusters

start_cell <- rownames(ms_df)[which.max(ms_df$MS_2)]
cat("Start cell:", start_cell, "\n")

pdf(file="./Palantir/1-FHF2_dataSub_Palantir_ms-dm_ManualrootCell.pdf", width=10, height=10)
ggplot(ms_df, aes(MS_1, MS_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_point(
    data  = ms_df[start_cell, , drop = FALSE],
    color = "red", size = 5, shape = 18
  ) +
  geom_label(
    data = ms_df[start_cell, , drop = FALSE],
    aes(x = MS_1, y = MS_2, label = "Start"),
    color = "red", vjust = -1, fontface = "bold",
    inherit.aes = FALSE
  ) +
  theme_classic()
dev.off()

# Identify tip cells manually — one per branch
tip1 <- rownames(ms_df)[which.max(ms_df$MS_1)]
tip2 <- rownames(ms_df)[which.min(ms_df$MS_2)]

cat("Tip 1 (MS_1 branch):", tip1, "\n")
cat("Tip 2 (MS_2 branch):", tip2, "\n")

# Visualise both tips before running
pdf(file="./Palantir/1-FHF2_dataSub_Palantir_tips.pdf", width=10, height=10)
ggplot(ms_df, aes(MS_1, MS_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_point(data = ms_df[c(start_cell, tip1, tip2), , drop = FALSE],
             aes(x = MS_1, y = MS_2),
             color = "red", size = 5, shape = 18, inherit.aes = FALSE) +
  geom_label(data = ms_df[c(start_cell, tip1, tip2), , drop = FALSE],
             aes(x = MS_1, y = MS_2,
                 label = c("Start", "Tip1", "Tip2")),
             color = "red", vjust = -1, fontface = "bold",
             inherit.aes = FALSE) +
  theme_classic()
dev.off()

# Palantir pseudotime
seurat_obj <- Palantir.Pseudotime(
  seurat_obj,
  start_cell = start_cell,
  terminal_states = c("fate1" = tip1, "fate2" = tip2)
)

ps <- seurat_obj@misc$Palantir$Pseudotime
print(colnames(ps))
print(head(ps))

# Extract pseudotime results
ps <- seurat_obj@misc$Palantir$Pseudotime
print(head(ps))

# Add to meta.data for plotting
seurat_obj@meta.data[, colnames(ps)] <- ps

pdf(file="./Palantir/1-FHF2_dataSub_Palantir_FateProbabilities_manualTips.pdf", width=10, height=10)
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "ms",
  cols      = list(Entropy = "D")   # "D" = diverging palette for entropy
)

# On original UMAP
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "umap",
  cols      = list(Entropy = "D")
)
dev.off()

seurat_obj@meta.data$Pseudotime_fate1 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$fate1
seurat_obj@meta.data$Pseudotime_fate2 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$fate2

pdf("./Palantir/1-FHF2_dataSub_Palantir_FatePseudotimePlot.pdf", width=10, height=10)
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate1‑weighted pseudotime"))
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate2‑weighted‑weighted pseudotime"))
dev.off()

saveRDS(seurat_obj, file = "./Palantir/1-FHF2_dataSub_Palantir.RDS")


	#######
	# JCF #
	#######
JCF1_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
JCF1_dataSub <- subset(x = JCF1_dataPreSub, subset = seurat_clusters %in% c(1,2,4,5,8,11), invert=FALSE)

seurat_obj <- Palantir.RunDM(
  JCF1_dataSub,
  reduction = "harmony",
  n_components = 30
)

# Inspect new reductions added
seurat_obj@reductions

# Visualise the multiscale space
pdf(file="./Palantir/2-JCF1_dataSub_Palantir_ms-dm.pdf", width=10, height=10)
DimPlot2(seurat_obj, reduction = "ms", group.by = "seurat_clusters", label = TRUE)
DimPlot2(seurat_obj, reduction = "dm", group.by = "seurat_clusters", label = TRUE)
dev.off()

ms_df <- as.data.frame(Embeddings(seurat_obj, "ms")[, 1:2])
colnames(ms_df) <- c("MS_1", "MS_2")
ms_df$cluster   <- seurat_obj$seurat_clusters

# Pick tip cell as root
start_cell <- rownames(ms_df)[which.min(ms_df$MS_2)]
cat("Start cell:", start_cell, "\n")

pdf(file="./Palantir/2-JCF1_dataSub_Palantir_ms-dm_ManualrootCell.pdf", width=10, height=10)
ggplot(ms_df, aes(MS_1, MS_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_point(
    data  = ms_df[start_cell, , drop = FALSE],
    color = "red", size = 5, shape = 18
  ) +
  geom_label(
    data = ms_df[start_cell, , drop = FALSE],
    aes(x = MS_1, y = MS_2, label = "Start"),
    color = "red", vjust = -1, fontface = "bold",
    inherit.aes = FALSE
  ) +
  theme_classic()
dev.off()

# Palantir pseudotime
seurat_obj <- Palantir.Pseudotime(
  seurat_obj,
  start_cell = start_cell
)

ps <- seurat_obj@misc$Palantir$Pseudotime
print(colnames(ps))
print(head(ps))

# Add to meta.data for plotting
seurat_obj@meta.data[, colnames(ps)] <- ps

# On ms space
pdf(file="./Palantir/2-JCF1_dataSub_Palantir_FateProbabilities.pdf", width=10, height=10)
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "ms",
  cols      = list(Entropy = "D")   # "D" = diverging palette for entropy
)

# On original UMAP
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "umap",
  cols      = list(Entropy = "D")
)
dev.off()

seurat_obj@meta.data$Pseudotime_fate1 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$cell_84601
seurat_obj@meta.data$Pseudotime_fate2 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$cell_12041

pdf("./Palantir/2-JCF1_dataSub_Palantir_FatePseudotimePlot.pdf", width=10, height=10)
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate1‑weighted pseudotime"))
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate2‑weighted‑weighted pseudotime"))
dev.off()

saveRDS(seurat_obj, file = "./Palantir/2-JCF1_dataSub_Palantir.RDS")


	#######
	# SHF #
	#######
SHF2_dataPreSub <- subset(x = Merge3_Krup_Tyser_Pijuan_harmony2_res1.0, subset = TimeStampMerged %in% c("E7.25","E7.5","E7.75","E8.0","E8.25","E8.5"), invert=FALSE)
SHF2_dataSub <- subset(x = SHF2_dataPreSub, subset = seurat_clusters %in% c(1,6,9), invert=FALSE)


seurat_obj <- Palantir.RunDM(
  SHF2_dataSub,
  reduction = "harmony",
  n_components = 30
)

# Inspect new reductions added
seurat_obj@reductions

# Visualise the multiscale space
pdf(file="./Palantir/3-SHF2_dataSub_Palantir_ms-dm.pdf", width=10, height=10)
DimPlot2(seurat_obj, reduction = "ms", group.by = "seurat_clusters", label = TRUE)
DimPlot2(seurat_obj, reduction = "dm", group.by = "seurat_clusters", label = TRUE)
dev.off()

ms_df <- as.data.frame(Embeddings(seurat_obj, "ms")[, 1:2])
colnames(ms_df) <- c("MS_1", "MS_2")
ms_df$cluster   <- seurat_obj$seurat_clusters

# Pick tip cell as root
start_cell <- rownames(ms_df)[which.min(ms_df$MS_1)]
cat("Start cell:", start_cell, "\n")

pdf(file="./Palantir/3-SHF2_dataSub_Palantir_ms-dm_ManualrootCell.pdf", width=10, height=10)
ggplot(ms_df, aes(MS_1, MS_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_point(
    data  = ms_df[start_cell, , drop = FALSE],
    color = "red", size = 5, shape = 18
  ) +
  geom_label(
    data = ms_df[start_cell, , drop = FALSE],
    aes(x = MS_1, y = MS_2, label = "Start"),
    color = "red", vjust = -1, fontface = "bold",
    inherit.aes = FALSE
  ) +
  theme_classic()
dev.off()

# Identify tip cells manually — one per branch
tip1 <- rownames(ms_df)[which.max(ms_df$MS_2)]
tip2 <- rownames(ms_df)[which.min(ms_df$MS_2)]

cat("Tip 1 (MS_1 branch):", tip1, "\n")
cat("Tip 2 (MS_2 branch):", tip2, "\n")

# Visualise both tips before running
pdf(file="./Palantir/3-SHF2_dataSub_Palantir_tips.pdf", width=10, height=10)
ggplot(ms_df, aes(MS_1, MS_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_point(data = ms_df[c(start_cell, tip1, tip2), , drop = FALSE],
             aes(x = MS_1, y = MS_2),
             color = "red", size = 5, shape = 18, inherit.aes = FALSE) +
  geom_label(data = ms_df[c(start_cell, tip1, tip2), , drop = FALSE],
             aes(x = MS_1, y = MS_2,
                 label = c("Start", "Tip1", "Tip2")),
             color = "red", vjust = -1, fontface = "bold",
             inherit.aes = FALSE) +
  theme_classic()
dev.off()

# Palantir pseudotime
seurat_obj <- Palantir.Pseudotime(
  seurat_obj,
  start_cell = start_cell,
  terminal_states = c("fate1" = tip1, "fate2" = tip2)
)

ps <- seurat_obj@misc$Palantir$Pseudotime
print(colnames(ps))
print(head(ps))

# Extract pseudotime results
ps <- seurat_obj@misc$Palantir$Pseudotime
print(head(ps))

# Add to meta.data for plotting
seurat_obj@meta.data[, colnames(ps)] <- ps

pdf(file="./Palantir/3-SHF2_dataSub_Palantir_FateProbabilities_manualTips.pdf", width=10, height=10)
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "ms",
  cols      = list(Entropy = "D")   # "D" = diverging palette for entropy
)

# On original UMAP
DimPlot2(
  seurat_obj,
  features  = colnames(ps),
  reduction = "umap",
  cols      = list(Entropy = "D")
)
dev.off()

seurat_obj@meta.data$Pseudotime_fate1 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$fate1
seurat_obj@meta.data$Pseudotime_fate2 <- seurat_obj@meta.data$Pseudotime * seurat_obj@meta.data$fate2

pdf("./Palantir/3-SHF2_dataSub_Palantir_FatePseudotimePlot.pdf", width=10, height=10)
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate1", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate1‑weighted pseudotime"))
print(FeaturePlot(seurat_obj, features = "Pseudotime_fate2", reduction = "umap", order = TRUE, pt.size = 0.8) + scale_color_gradientn(colors = c("#0D0887FF","#CC4678FF","#F0F921FF"), na.value = "#E4E4E4") + labs(x="UMAP1", y="UMAP2", title= "Fate2‑weighted‑weighted pseudotime"))
dev.off()

saveRDS(seurat_obj, file = "./Palantir/3-SHF2_dataSub_Palantir.RDS")


#=======================================================================================================================================================================#
###################################################========      TSCAN & Palantir pseudotimes comparison vs. Slingshot      ========#####################################
#=======================================================================================================================================================================#
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
  library(patchwork)

	# Load data
FHF2_dataSub_TSCAN <- readRDS(file = "./TSCAN/1-FHF2_dataSub_TSCAN.RDS")
JCF1_dataSub_TSCAN <- readRDS(file = "./TSCAN/2-JCF1_dataSub_TSCAN.RDS")
SHF2_dataSub_TSCAN <- readRDS(file = "./TSCAN/3-SHF2_dataSub_TSCAN.RDS")

FHF2_dataSub_Palantir <- readRDS(file = "./Palantir/1-FHF2_dataSub_Palantir.RDS")
JCF1_dataSub_Palantir <- readRDS(file = "./Palantir/2-JCF1_dataSub_Palantir.RDS")
SHF2_dataSub_Palantir <- readRDS(file = "./Palantir/3-SHF2_dataSub_Palantir.RDS")

# Helper function
compare_pseudotime <- function(ref_col, tool_obj, tool_col, lineage_name, tool_name, color = "steelblue") {
  
  # Extract reference pseudotime, remove NAs
  ref_pt <- data.frame(
    Pseudotime_ref = dataset@meta.data[[ref_col]],
    row.names = rownames(dataset@meta.data)
  )
  ref_pt <- ref_pt[!is.na(ref_pt$Pseudotime_ref), , drop = FALSE]
  
  # Extract tool pseudotime using specified column
  pt_tool <- data.frame(
    Pseudotime_tool = tool_obj@meta.data[[tool_col]],
    row.names = rownames(tool_obj@meta.data)
  )
  
  # Keep only cells present in ref_pt before merging
  pt_tool <- pt_tool[rownames(pt_tool) %in% rownames(ref_pt), , drop = FALSE]
  
  # Merge on common cells
  merged <- merge(ref_pt, pt_tool, by = "row.names")
  rownames(merged) <- merged$Row.names; merged$Row.names <- NULL
  
  cat("\n──", lineage_name, "-", tool_name, "──\n")
  cat("Cells in Slingshot reference:", nrow(ref_pt), "\n")
  cat("Cells after filtering:", nrow(merged), "\n")
  
  # Pearson correlation
  cor_res <- cor.test(merged$Pseudotime_ref, merged$Pseudotime_tool, method = "pearson")
  cat("r =", round(cor_res$estimate, 3), "| p =", formatC(cor_res$p.value, format = "e", digits = 2), "\n")
  
  # Scatter plot
  p <- ggplot(merged, aes(x = Pseudotime_ref, y = Pseudotime_tool)) +
    geom_point(alpha = 0.3, size = 0.8, color = color) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(cor_res$estimate, 3),
                            "\np = ", formatC(cor_res$p.value, format = "e", digits = 2))) +
    labs(x = "Slingshot pseudotime",
         y = paste(tool_name, "pseudotime"),
         title = paste0(lineage_name, " — ", tool_name)) +
    theme_classic()
  
  list(plot = p, merged = merged, cor = cor_res)
}

# ── Run for each lineage x tool ───────────────────────────────────────────────
res_FHF_TSCAN     <- compare_pseudotime("FHF_Pseudotime",  FHF2_dataSub_TSCAN,    "TSCAN_path1", "FHF",  "TSCAN",    color = "steelblue")
res_FHF_Palantir  <- compare_pseudotime("FHF_Pseudotime",  FHF2_dataSub_Palantir, "Pseudotime_fate1", "FHF",  "Palantir", color = "darkorange")
res_JCF_TSCAN     <- compare_pseudotime("JCF_Pseudotime",  JCF1_dataSub_TSCAN,    "TSCAN_path2", "JCF",  "TSCAN",    color = "steelblue")
res_JCF_Palantir  <- compare_pseudotime("JCF_Pseudotime",  JCF1_dataSub_Palantir, "Pseudotime_fate2", "JCF",  "Palantir", color = "darkorange")
res_aSHF_TSCAN    <- compare_pseudotime("aSHF_Pseudotime", SHF2_dataSub_TSCAN,    "TSCAN_path1", "aSHF", "TSCAN",    color = "steelblue")
res_aSHF_Palantir <- compare_pseudotime("aSHF_Pseudotime", SHF2_dataSub_Palantir, "Pseudotime_fate2", "aSHF", "Palantir", color = "darkorange")
res_pSHF_TSCAN    <- compare_pseudotime("pSHF_Pseudotime", SHF2_dataSub_TSCAN,    "TSCAN_path2", "pSHF", "TSCAN",    color = "steelblue")
res_pSHF_Palantir <- compare_pseudotime("pSHF_Pseudotime", SHF2_dataSub_Palantir, "Pseudotime_fate1", "pSHF", "Palantir", color = "darkorange")

# ── Save to PDF ───────────────────────────────────────────────────────────────
pdf("./3-Pseudotime_comparison.pdf",
    width = 12, height = 5)

(res_FHF_TSCAN$plot  | res_FHF_Palantir$plot)  + plot_annotation(title = "FHF")
(res_JCF_TSCAN$plot  | res_JCF_Palantir$plot)  + plot_annotation(title = "JCF")
(res_aSHF_TSCAN$plot | res_aSHF_Palantir$plot) + plot_annotation(title = "aSHF")
(res_pSHF_TSCAN$plot | res_pSHF_Palantir$plot) + plot_annotation(title = "pSHF")

dev.off()

pdf("./3-Pseudotime_comparison2.v2.pdf",
    width = 20, height = 10)

(res_FHF_TSCAN$plot | res_JCF_TSCAN$plot | res_aSHF_TSCAN$plot | res_pSHF_TSCAN$plot) /
(res_FHF_Palantir$plot | res_JCF_Palantir$plot | res_aSHF_Palantir$plot | res_pSHF_Palantir$plot)

dev.off()
