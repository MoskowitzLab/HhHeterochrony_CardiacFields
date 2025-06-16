AddDoublets <- function(seuratObj, scrubFile){
  
  scrub <- as.data.frame(read.table(scrubFile, sep = '\t', header = TRUE))
  scrubSet <- scrub[scrub$cellID %in% colnames(seuratObj), ]
  
  seuratObj <- AddMetaData(
    object = seuratObj,
    metadata = scrubSet$doublet_score,
    col.name = 'doublet_score'
  )
  
  seuratObj <- AddMetaData(
    object = seuratObj,
    metadata = scrubSet$predicted_doublet,
    col.name = 'doublet_calls'
  )

  return(seuratObj)
  
}
