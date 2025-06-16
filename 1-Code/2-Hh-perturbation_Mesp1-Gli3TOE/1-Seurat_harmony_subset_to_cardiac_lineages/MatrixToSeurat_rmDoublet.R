
	
MatrixToSeurat_rmDoublet <- function(DGEmatrix, whatSample, whatGenotype, doubletfile){
	
	  seuratObject <- CreateSeuratObject(
	  counts = DGEmatrix,
	  project = "HH-ScRNAseq",
	  assay = "RNA",
	  min.cells = 5,
	  min.features = 500,
	)
	
	seuratObject$Sample <- whatSample
	seuratObject$Genotype <- whatGenotype
	
	seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^mt-")
	Idents(seuratObject) <- seuratObject$Sample
	
	source("/project/imoskowitz/shared/data_collaborators/AlexG/code/AddDoublets.R")
	seuratObject <- AddDoublets(seuratObject, doubletfile) 
	
	print("Number of doublets detected")
	print(length(subset(seuratObject$doublet_calls, seuratObject$doublet_calls == "True")))
	
	print("Percent of cells called as doublets")
	print(length(subset(seuratObject$doublet_calls, seuratObject$doublet_calls == "True"))/length(seuratObject$doublet_calls)*100)
	
	
	
	seuratObject <- subset(seuratObject,
	              subset = doublet_calls == "False")
	
	print("Percent of cells passing QC metrics")
	print(ncol(seuratObject))


	return(seuratObject)

}
