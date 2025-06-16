matchIDToSymbol2 <- function(DGEmatrix, genenames){
  
  #name columns in geneID to Symbol conversion file
  colnames(genenames) <- c("GeneID", "Type", "Symbol")  
  
  #add the Symbol column - of the intial dge only first 2 columns are used eles it's too slow for last dfs
  added_annot_DF <- DGEmatrix[, 1:2] %>%  
    inner_join(genenames, ., by= c("GeneID" =  "gene"))
  
  #the Df contains redundant entries for Symbol. These are removed by choosing the first entry
  #for any redundant symbol
  added_annot_DF <- added_annot_DF %>%
    group_by(Symbol) %>%
    summarise_all(dplyr::first) %>%
    as.data.frame()
  
  
  added_annot_DF <- DGEmatrix %>%  
    inner_join(added_annot_DF[,1:3], ., by= c("GeneID" =  "gene"))
  
  #add symbol as rownames
  rownames(added_annot_DF) <- added_annot_DF$Symbol
  
  #remove annotation columns (1 and 2)
  added_annot_DF <- added_annot_DF[ , -(1:3)]
  
  #return dataframe with 
  return(added_annot_DF)
  
}


