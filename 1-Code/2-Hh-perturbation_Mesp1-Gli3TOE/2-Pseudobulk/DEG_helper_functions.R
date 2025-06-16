# TXI import # this is actually Carlos' function
txi_fromtximport_salmon <- function(salmonFolder, 
                                    metaTable) {
  
  require(DESeq2)
  require(tximport)
  require(dplyr)
  require(stringr)
  require(readr)
  
  star_salmon_path = salmonFolder
  tx2gene = read.delim( paste0(star_salmon_path, "/tx2gene.tsv"))
  metad = metaTable
  rownames(metad) <- metad$Sample
  files <- file.path(star_salmon_path, metad$Sample, "quant.sf")
  toremove <- c(star_salmon_path,"/", "quant.sf")
  names(files) <- str_remove_all(files, paste(toremove, collapse = "|"))
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  return(txi)
}

################### RESULTS TRANSFORMATION ########################
results_transform_joinCounts <- function(res, normCounts, by="gene_id", alpha=0.05, lfcThreshold=0){
  res <- res %>%
    as.data.frame() %>%
    rownames_to_column(by) %>%
    .[complete.cases(.$padj),] %>%
    mutate(
      significance = ifelse(padj < alpha, TRUE, FALSE),
      differentialExpression = case_when(padj < alpha  & log2FoldChange > lfcThreshold ~ "up",
                                         padj < alpha & log2FoldChange < -lfcThreshold ~ "down",
                                         TRUE ~ "not significant"
      ) %>%
    left_join(normCounts, by = by)
    )
  return(res)
}

results_transform <- function(res, alpha=0.05, lfcThreshold=0, rownames_to_column = "gene_id"){

    res <- res %>%
    as.data.frame() %>%
    rownames_to_column(rownames_to_column) %>%
    .[complete.cases(.$padj),] %>%
    mutate(
      significance = ifelse(padj < alpha & abs(log2FoldChange) > lfcThreshold, TRUE, FALSE),
      differentialExpression = case_when(padj < alpha  & log2FoldChange > lfcThreshold ~ "up",
                                         padj < alpha & log2FoldChange < -lfcThreshold ~ "down",
                                         TRUE ~ "not significant"
      )
    )
  return(res)
}


####### VOLCANO PLOT #######
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
plotVolcano <- function(results,
                        p_val = "padj",
                        title = NULL,
                        color_column = "differentialExpression",
                        annotations_column = "label",
                        point_colors = c("steelblue", "gray88", "coral3"),
                        lfcThreshold = 0,
                        padjThreshold = 0.05,
                        use_asinh_scale_x = FALSE,
                        use_asinh_scale_y = FALSE,
                        add_threshold_lines = TRUE,
                        show_degs_number = TRUE,
                        annotation_coords = c(-max_lfc / 1.1, max_lfc /
                                                1.1, max_padj/2),
                        show_labels = TRUE)  {
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  require(ggrastr)
  
  max_padj <-
    range(-log10(results$padj), finite = TRUE, na.rm = TRUE)[2]
  padj_breaks <-
    seq(0, max_padj, length.out = 5) %>% round_any(accuracy = 10, x = .)
  
  max_lfc <-
    range(results$log2FoldChange,
          finite = TRUE,
          na.rm = TRUE) %>% max()
  lfc_breaks <-
    seq(-max_lfc, max_lfc, length.out = 8) %>% round_any(accuracy = 2, x = .)
  
  ups <- sum(results$differentialExpression == "up")
  downs <- sum(results$differentialExpression == "down")
  
  gg <-
    ggplot(
      results,
      aes(
        x = log2FoldChange,
        y = -log10(!!sym(p_val)),
        col = !!sym(color_column),
        label = !!sym(annotations_column)
      )
    ) +
    ggrastr::rasterise(geom_point(size = 1.2, alpha = 0.9), dpi=300) +
    
    scale_color_manual(
      values = point_colors,
      labels = c("Downregulated", "Not signicant", "Upregulated")
    ) +
    
    theme_classic() +
    theme(
      text = element_text(family = "sans"),
      legend.position = 'bottom',
      legend.direction = 'horizontal',
      legend.title = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    
    labs(x = expression(paste("log"[2] , " Fold Change")),
         y = expression(paste("-log"[10], " Adj. Pval")))
  
  if (use_asinh_scale_y) {
    gg <-
      gg + scale_y_continuous(
        trans = 'asinh',
        limits = c(0, max_padj),
        breaks = c(padj_breaks)
      )
  } else {
    gg <- gg + scale_y_continuous()
  }
  
  if (use_asinh_scale_x) {
    gg <-
      gg + scale_x_continuous(
        trans = "asinh",
        limits = c(-max_lfc, max_lfc),
        breaks = c(lfc_breaks)
      )
  } else {
    gg <- gg + scale_x_continuous()
  }
  
  if (add_threshold_lines) {
    gg <- gg +
      geom_hline(
        yintercept = -log10(padjThreshold),
        linetype = "dashed",
        color = "gray61"
      ) +
      geom_vline(
        xintercept = c(lfcThreshold,-lfcThreshold),
        linetype = 'dashed',
        color = 'gray61'
      )
  }
  
  if (show_labels) {
    gg <- gg +
      geom_label_repel(
        box.padding = 0.25,
        point.padding = 0.5,
        force_pull = 20,
        size = 3,
        show.legend = FALSE,
        max.overlaps = Inf
      )
    
    if (show_degs_number) {
      gg <- gg +
        annotate(
          "text",
          x = c(annotation_coords[1], annotation_coords[2]),
          y = annotation_coords[3],
          color = c(point_colors[1], point_colors[3]),
          label = c(downs, ups),
          size = c(5, 5),
          fontface = c("bold", "bold")
        )
    }
  }
  
  gg <- gg + ggtitle(title)
  
  return(gg)
}
label_top_degs <- function(df, genes = "gene_id", number = 20) {
  topUP <- df %>% 
    subset(differentialExpression == "up") %>% 
    slice_min(n = number, order_by = padj) %>%
    head(n=number)
  topDOWN <- df %>% 
    subset(differentialExpression == "down") %>% 
    slice_min(n = number, order_by = padj) %>%
    head(n=number)
  
  df$label <- ifelse(df[[genes]] %in% c(topUP[[genes]], topDOWN[[genes]]), df[[genes]], NA)
  return(df)
}


label_top_degs_lfc <- function(df, genes = "gene_id", number = 20) {
  topUP <- df %>% 
    subset(differentialExpression == "up") %>% 
    slice_max(n = number, order_by = log2FoldChange) %>%
    head(n = number)
  topDOWN <- df %>% 
    subset(differentialExpression == "down") %>% 
    slice_max(n = number, order_by = -log2FoldChange) %>%
    head(n = number)
  
  df$label <- ifelse(df[[genes]] %in% c(topUP[[genes]], topDOWN[[genes]]), df[[genes]], NA)
  return(df)
}


top_degs_annotation_merged <- function(df, genes = "gene_id", number = 40, suffix = c(".x", ".y"), annotateLFC = TRUE) {
  x <- substr(suffix[1], 2, nchar(suffix[1]))
  y <- substr(suffix[2], 2, nchar(suffix[2]))
  padj.x <- paste0("padj", suffix[1])
  padj.y <- paste0("padj", suffix[2])
  log2FoldChange.x <- paste0("log2FoldChange", suffix[1])
  log2FoldChange.y <- paste0("log2FoldChange", suffix[2])
  
  df_annotated <- df %>%
    group_by(Significance) %>%
    mutate(
      annotations = case_when(
        Significance == "Significant in both" ~ .data[[genes]],
        Significance == paste0("Significant in ", x) ~ .data[[genes]],
        Significance == paste0("Significant in ", y) ~ .data[[genes]],
        TRUE ~ NA_character_
      )
    )
  
  
  genes_to_annotate <- df_annotated %>%
    filter(!is.na(annotations)) %>%
    group_by(Significance) %>%
    arrange(across(all_of(paste0("padj", suffix)))) %>%
    slice_min(n = as.integer(number / 3), order_by = across(all_of(paste0("padj", suffix)))) %>%
    pull(annotations)
  
  if (annotateLFC) {
    high_LFC_genes <- df_annotated %>%
      filter(abs(!!sym(log2FoldChange.x)) > 10 | abs(!!sym(log2FoldChange.y)) > 10) %>%
      pull({{genes}})
    
    genes_to_annotate <- union(genes_to_annotate, high_LFC_genes)
  }

  
  df$annotations <- ifelse(df[[genes]] %in% genes_to_annotate, df[[genes]], NA)
  return(df)
}





transformMerged <- function(merged_df, suffix = c(".x", ".y")){
  x <- substr(suffix[1], 2, nchar(suffix[1]))
  y <- substr(suffix[2], 2, nchar(suffix[2]))
  
  significance.x <- paste0("significance", suffix[1])
  significance.y <- paste0("significance", suffix[2])
  Log2FoldChange.x <- paste0("log2FoldChange", suffix[1])
  Log2FoldChange.y <- paste0("log2FoldChange", suffix[2])
                                     
  merged_df <- merged_df %>%
    mutate(
      Significance = case_when(
        .data[[significance.x]] == TRUE & .data[[significance.y]] == TRUE ~ "Significant in both",
        .data[[significance.x]] == TRUE & .data[[significance.y]] == FALSE ~ paste0("Significant in ", x),
        .data[[significance.x]] == FALSE & .data[[significance.y]] == TRUE ~ paste0("Significant in ", y),
        TRUE ~ "Not significant in either"
      ),
      Quadrant = case_when(
        .data[[Log2FoldChange.x]] > 0 & .data[[Log2FoldChange.y]] >=0 ~ "upper-right",
        .data[[Log2FoldChange.x]] <0 & .data[[Log2FoldChange.y]] >=0 ~ "upper-left",
        .data[[Log2FoldChange.x]] <0  & .data[[Log2FoldChange.y]] <=0  ~ "lower-left",
        .data[[Log2FoldChange.x]] >0 & .data[[Log2FoldChange.y]] <=0  ~ "lower-right"
      )
    )
  return(merged_df)
}


######## GENE ONTOLOGY ################
runGO <- function(test.vector, 
                  univr_vector, 
                  ont="BP",
                  organism='mouse',
                  simplify=TRUE,
                  keyType = "ENTREZID") {
  
  require(dplyr)
  require(ggplot2)
  require(clusterProfiler)
  require(org.Mm.eg.db)
  #require(org.Hs.eg.db)
  
  if(organism=='mouse') {
    OrgDb= org.Mm.eg.db
  } else {
    OrgDb=org.Hs.eg.db
  }
  
  
  goPATH <- enrichGO(gene=test.vector, 
                     OrgDb = OrgDb, 
                     universe=univr_vector,
                     ont=ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable = TRUE,
                     keyType = keyType)
  
  #goPATH <- setReadable(goPATH, OrgDb = OrgDb)  
  
  if (simplify) {
    goPATH <- simplify(goPATH)
  }
  
  goPATH@result <-  goPATH@geneSets %>% 
    lapply(. ,length) %>% 
    unlist %>% data.frame("GenesInGo" = .) %>%
    tibble::rownames_to_column("ID") %>% 
    left_join(goPATH@result, .)
  
  return(goPATH)
  
}

barplotGO <- function(object, topTerms = 10, barColor = "gray50", outlineColor = NULL, subtitle = NULL, strWrap = 40, barWidth = 0.8, yPosition = "left", alpha = 1) {
  require(ggplot2)
  require(dplyr)
  require(scales)
  goPATH <- data.frame(object) %>%
    arrange(p.adjust) %>%
    head(topTerms)
  
  max_padj <- -log10(goPATH$p.adjust) %>% max()
  
  gg <- ggplot(data = goPATH, aes(x = -log10(p.adjust), y = Description)) +
    geom_bar(stat = "identity", fill = barColor, width = barWidth, alpha = alpha) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = 'transparent'),
      plot.subtitle = element_text(size = 20, face = "bold", color = barColor, family = "sans", hjust = ifelse(yPosition == "left", 1, 0)),
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = .3), 
      axis.ticks.x = element_blank(),
      text = element_text(family = "sans", color = outlineColor),
      axis.title.x = element_text(size = 20, hjust = 0.5),
      axis.text.y = element_text(size = 20),
      axis.text.x = element_text(size = 18),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      x = expression(paste("-log"[10], " Adj. Pval.")),
      y = "",
      subtitle = subtitle
    ) +
    scale_y_discrete(
      position = yPosition, 
      labels = wrap_format(strWrap),
      limits = as.vector(goPATH[order(-log10(goPATH$p.adjust), decreasing = FALSE), ])$Description, 
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, max_padj))
  
  if (yPosition == "left") {
    gg <- gg + 
      scale_x_reverse(expand = expansion(mult = c(0, 0)), limits = c(max_padj, 0))
  }
  
  return(gg)
}



plotCorrelations <- function(data, 
                             suffix = c(".x", ".y"),
                             title = NULL, 
                             color_column = "Significance", 
                             annotations_column = "annotations", 
                             point_colors = c("gray88",
                                              "darkseagreen",
                                              "plum2",
                                              "lightsalmon"),
                             show_labels = TRUE) {
  
  require(ggplot2)
  require(ggrepel)
  require(ggpmisc)
  require(dplyr)
  
  # Create summary data frame for annotations
  summary_data <- data %>%
    group_by(Significance, Quadrant) %>%
    summarise(count = n(), .groups = 'drop') %>%
    filter(Significance != "Not significant in either") # Filter out "Not significant in either"
  
  # Create coordinates for placing the text
  summary_data <- summary_data %>%
    mutate(x_coord = case_when(
      Quadrant == "upper-right" ~ max(data[[paste0("log2FoldChange", suffix[1])]]) * 0.75,
      Quadrant == "upper-left" ~ min(data[[paste0("log2FoldChange", suffix[1])]]) * 0.75,
      Quadrant == "lower-right" ~ max(data[[paste0("log2FoldChange", suffix[1])]]) * 0.75,
      Quadrant == "lower-left" ~ min(data[[paste0("log2FoldChange", suffix[1])]]) * 0.75,
      TRUE ~ 0
    ),
    y_coord = case_when(
      Quadrant == "upper-right" ~ max(data[[paste0("log2FoldChange", suffix[2])]]) * 0.65, # Adjust y_coord lower
      Quadrant == "upper-left" ~ max(data[[paste0("log2FoldChange", suffix[2])]]) * 0.65,  # Adjust y_coord lower
      Quadrant == "lower-right" ~ min(data[[paste0("log2FoldChange", suffix[2])]]) * 0.65, # Adjust y_coord lower
      Quadrant == "lower-left" ~ min(data[[paste0("log2FoldChange", suffix[2])]]) * 0.65,  # Adjust y_coord lower
      TRUE ~ 0
    )) %>%
    mutate(offset = case_when(
      Significance == "Significant in both" ~ 1,
      Significance == "Significant in DoxD18/CtlD18" ~ 2,
      Significance == "Significant in DoxD25/CtlD25" ~ 3,
      TRUE ~ 0
    )) %>%
    mutate(y_coord = y_coord + offset * 0.05 * diff(range(data[[paste0("log2FoldChange", suffix[2])]])))
  
  coords <- c(min(data[[paste0("log2FoldChange", suffix[1])]]) + 1, max(data[[paste0("log2FoldChange", suffix[2])]]))
  correlation <- cor(data[[paste0("log2FoldChange", suffix[1])]], data[[paste0("log2FoldChange", suffix[2])]], use = "pairwise.complete.obs")
  
  gg <- ggplot(data, aes(x = !!sym(paste0("log2FoldChange", suffix[1])), y = !!sym(paste0("log2FoldChange", suffix[2])), col = !!sym(color_column), label = !!sym(annotations_column))) +
    geom_point(size=1, alpha = 0.8) +
    ggtitle(title) +
    scale_color_manual(values = point_colors) +
    theme_classic() +
    theme(text = element_text(family = "sans"),
          legend.position = 'right', legend.direction = 'vertical', 
          axis.title = element_text(size=12),
          plot.title = element_text(size=13, face = "bold")) +
    labs(title = title,
         x= paste0("LFC", suffix[1]),
         y= paste0("LFC", suffix[2])) +
    geom_quadrant_lines()
  

  if (show_labels) {
    gg <- gg +
      geom_label_repel(
        box.padding = 0.25,
        point.padding = 0.5,
        force_pull = 20,
        size = 3,
        show.legend = FALSE,
        max.overlaps = Inf)
  }
  
  gg <- gg +
    geom_text(data = summary_data, aes(x = x_coord, y = y_coord, label = count, color = Significance), size = 4, fontface = "bold", show.legend = FALSE)
  
  return(gg)
}


vennDiagram <- function(x=list(), category.names, fill = c("#B6D8F2", "#F7CFDF"),main=NULL, savefile = FALSE, imagetype = "svg", filename = NULL, ext.text=TRUE, cat.pos = c(-50,50)) {
  require(VennDiagram)
  require(svglite)
  grid.newpage()
  venn_plot <- venn.diagram(
    x = x,
    filename=NULL,
    category.names = category.names,
    fill = fill,
    main = main,
    main.fontfamily	= "Arial Narrow",
    main.cex = 2.2,
    main.pos = c(0.5, .85),
    alpha = 0.6,
    lwd = 0,
    cat.fontcolor = c("black", "black"),
    cat.cex = 2,
    cat.fontfamily = "Arial Narrow",
    cex = 2,
    fontfamily = "Arial Narrow",
    cat.dist = 0.05,
    margin = 0.1,
    label.col = "black",
    label.cex = 1.5,
    category.percentage = TRUE,
    scaled = TRUE,
    print.mode = c("raw", "percent"),
    euler.d = TRUE,
    units = "in",
    inverted = FALSE,
    ext.text = ext.text,
    cat.pos = cat.pos
  )
  
  if (savefile) {
    if (!is.null(filename)) {
      if (imagetype == "svg") {
        svglite(filename, width = 6, height = 6, bg="transparent")
        grid.draw(venn_plot)
        dev.off()
      } else if (imagetype == "pdf") {
        pdf(filename, width = 6, height = 6, bg="transparent")
        grid.draw(venn_plot)
        dev.off()
      } else {
        cat("Invalid imagetype. Supported types are 'svg' and 'pdf'.\n")
      }
    } else {
      cat("You must provide a valid filename when saving the diagram.\n")
    }
  } else {
    grid.draw(venn_plot)
  }
}

##### FILTER COUNTS #########
filter_sigCounts <- function(dds, metadata,  metadata_column, conditions, results) {
  filtered_dds <- dds[, dds[[metadata_column]] %in% conditions]
  filtered_metadata <- metadata[dds[[metadata_column]] %in% conditions, ]
  normCounts <- counts(filtered_dds, normalized=TRUE) %>% as.data.frame()
  sigNormCounts <- normCounts[subset(results, significance==TRUE)$gene_id, ]
  return(list(sigNormCounts = sigNormCounts, filtered_metadata = filtered_metadata))
}

filterCounts <- function(dds, metadata,  metadata_column, conditions, results) {
  filtered_dds <- dds[, dds[[metadata_column]] %in% conditions]
  filtered_metadata <- metadata[dds[[metadata_column]] %in% conditions, ]
  normCounts <- counts(filtered_dds, normalized=TRUE) %>% as.data.frame()
  return(list(normCounts = normCounts, filtered_metadata = filtered_metadata))
}
##### FILTER VSD #########
filter_vsd <- function(vsd, metadata, metadata_column, conditions) {
  filtered_vsd <- vsd[, vsd[[metadata_column]] %in% conditions]
  filtered_metadata <- metadata[vsd[[metadata_column]] %in% conditions, ]
  return(list(vsd = filtered_vsd, metadata = filtered_metadata))
}


############ BOX PLOTS ##############
# Function to create the plot for a specific gene and type
genes_boxplot <- function(data, x, y, perform_test = TRUE, test_method = "t.test", method_args = list(var.equal = TRUE), perform_test_manually = FALSE, stat_test_df = NULL, stat_test_label = NULL, colors = c(Ctl = "#B6D8F2", OE = "#F7CFDF"), yPosition = "left", subtitle = "") {
  
  require(ggplot2)
  require(dplyr)
  require(ggpubr)

  plot <- ggplot(data, aes(x = !!sym(x), y = !!sym(y))) +
    geom_boxplot(outlier.fill = NA, alpha = 0.9, aes(fill = !!sym(x))) +
    geom_jitter(aes(color = !!sym(x), fill = !!sym(x)), shape = 21, color = "black") +
    scale_fill_manual(values = colors) 
  
  if (perform_test) {
    plot <- plot + stat_compare_means(
      method = test_method,
      method.args = method_args,
      aes(group = !!sym(x)),
      label = "p.format",
      label.x.npc = "center",
      label.y.npc = "middle",
      family = "Arial Narrow",
      size = 3,
      geom = "label"
    ) 
  }
  
  if (perform_test_manually){
   plot <- plot + 
      stat_pvalue_manual(
        data = stat_test_df,
        label = stat_test_label,
        size=3,
        remove.bracket = F,
        family = "Arial Narrow",
        linetype=0
        )
    }
  
  plot <- plot + facet_wrap(~ .data$gene_id, scales = "free_y", ncol = 1, shrink = TRUE) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(family = "Arial Narrow"),
          legend.position = "none",
          panel.grid = element_blank())
  
  if (yPosition == "left") {
    plot <- plot + labs(subtitle = subtitle, x = "", y = expression(paste("log"[2], "TPM"))) + scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))
  } else {
    plot <- plot + labs(subtitle = subtitle, x = "", y = "")
  }
  
  if (yPosition == "right") {
    plot <- plot + scale_y_continuous(position = yPosition, expand = expansion(mult = c(0.2, 0.2))) 
  }
  
  return(plot)
}

### Extract genomic region from Homer's annotatePeaks.pl
extract_genomic_region <- function(annotation, patterns = c("promoter-TSS", "TTS", "Intergenic", "intron", "exon", "3' UTR", "5' UTR", "non-coding")) {
  
  # Loop through patterns and find matching one
  for (pattern in patterns) {
    if (startsWith(annotation, pattern)) {
      return(pattern) # Return the matched pattern
    }
  }
  
  return(NA) # Return NA if no match found
  
}
  
  
  
  
  
############ PAIRWISE KS TEST ##############
  pairwise_ks_test <- function(data, variable_column, group_column, p.adjust.method = "fdr") {
    # Define the groups
    groups <- unique(data[[group_column]])
    
    # Create an empty dataframe to store the results
    result_df <- data.frame(group1 = character(), group2 = character(), n1 = numeric(),
                            n2 = numeric(), p = numeric(), p.adj = numeric(), padj.signif = character(),
                            max_value_group1 = numeric(), max_value_group2 = numeric(),
                            stringsAsFactors = FALSE)
    
    # Perform pairwise Kolmogorov-Smirnov test and adjust p-values
    pvalues <- c()
    index <- 1
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        group1 <- groups[i]
        group2 <- groups[j]
        
        # Subset data for the two groups
        data_group1 <- data[data[[group_column]] == group1,]
        data_group2 <- data[data[[group_column]] == group2,]
        
        # Calculate the KS statistic
        ks_stat <- ks.test(data_group1[[variable_column]], data_group2[[variable_column]])
        pvalues[index] <- ks_stat$p.value
        
        # Number of observations in each group
        n1 <- length(data_group1[[variable_column]])
        n2 <- length(data_group2[[variable_column]])
        
        # Max value for each group
        max_value_group1 <- mean(data_group1[[variable_column]])#as.numeric(names(sort(table(data_group1[[variable_column]]), decreasing = TRUE))[1])
        max_value_group2 <- mean(data_group2[[variable_column]])#as.numeric(names(sort(table(data_group2[[variable_column]]), decreasing = TRUE))[1])
        
        # Store the results
        result_df <- rbind(result_df, data.frame(group1 = group1, group2 = group2,
                                                 n1 = n1, n2 = n2,
                                                 p = ks_stat$p.value, p.adj = 0,
                                                 padj.signif = 0,
                                                 max_value_group1 = max_value_group1,
                                                 max_value_group2 = max_value_group2))
        index <- index + 1
      }
    }
    padj <- p.adjust(pvalues, method = p.adjust.method)
    result_df$p.adj <- padj
    cutpoints <- c(0, 1e-04, 0.001, 0.01, 0.05, 1)
    symbols <- c("****", "***", "**", "*", "ns")
    result_df$padj.signif <- cut(result_df$p.adj, cutpoints, labels = symbols, include.lowest = TRUE)
    
    return(result_df)
  }

##### REMOVE LOWER VARIANCE
removeVar <- function(vsd, var_to_remove) {
  variable_genes <- assay(vsd) %>%
    rowVars() %>%
    order(decreasing=T) %>%
    head(n = length(vsd) * (1 - var_to_remove)) %>%
    vsd[. ,]
  return(variable_genes)
}
