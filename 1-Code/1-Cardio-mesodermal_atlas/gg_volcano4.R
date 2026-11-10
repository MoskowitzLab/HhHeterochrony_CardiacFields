# gg_volcano: Plot a volcano plot from a differential expression results dataframe.
gg_volcano4 <- function(df,
                       abs_lfc = 0,
                       p.cutoff = 0.05,
                       up.color = "steelblue3",
                       down.color = "steelblue3",
                       unSig.color = "grey88",
                       dotSize = 0.5,
                       alpha = 1,
                       theme_style = "minimal",
                       scale_axes = TRUE,
                       xlim = NULL,
                       ylim = NULL,
                       count_position = "bottom",
                       use_asinh = FALSE,
                       n_labels_up = 10,
                       n_labels_down = 10,
                       label_TF_only = FALSE,
                       tf_table_path = "/project/imoskowitz/shared/annotations/Mus_musculus_TF.txt",
                       align_labels = FALSE
) {
  require(dplyr)
  require(ggplot2)
  require(ggrastr)
  require(ggrepel)
  require(scales)

  # Detect log fold change and p-value columns
  lfc <- grep("log2FoldChange|logFC|lfc_|_lfc|avg_log2FC", colnames(df), ignore.case = T)
  p.col <- grep("padj|FDR|p_val_adj", colnames(df))

  if (length(lfc) == 0 || length(p.col) == 0) {
    stop("Could not find log fold change or adjusted p-value columns in the data frame")
  }

  # Prepare data
  df$logp <- -log10(df[[p.col]] + .Machine$double.xmin)
  df$lfc <- df[[lfc]]

  # Remove infinite values
  df <- df[is.finite(df$logp) & is.finite(df$lfc), ]

  # Calculate axis limits
  if (scale_axes) {
    maxlfc <- max(abs(range(df$lfc, finite = T, na.rm = T)))
    maxpadj <- max(df$logp, finite = T, na.rm = T)
    lfc_breaks <- pretty(c(-maxlfc, maxlfc), n = 6)
    padj_breaks <- pretty(c(0, maxpadj), n = 5)
    x_limits <- c(-maxlfc * 1.05, maxlfc * 1.05)
    y_limits <- c(0, maxpadj * 1.05)
  } else {
    if (is.null(xlim)) {
      maxlfc <- max(abs(range(df$lfc, finite = T, na.rm = T)))
      x_limits <- c(-maxlfc * 1.05, maxlfc * 1.05)
    } else {
      x_limits <- xlim
    }
    if (is.null(ylim)) {
      maxpadj <- max(df$logp, finite = T, na.rm = T)
      y_limits <- c(0, maxpadj * 1.05)
    } else {
      y_limits <- ylim
    }
    lfc_breaks <- pretty(x_limits, n = 6)
    padj_breaks <- pretty(y_limits, n = 5)
  }

  # Classify points
  Sig <- subset(df, df[[p.col]] < p.cutoff & abs(lfc) > abs_lfc)
  unSig <- subset(df, !rownames(df) %in% rownames(Sig))
  u <- subset(Sig, lfc > 0)
  d <- subset(Sig, lfc < 0)

  # Select theme
  theme_func <- switch(theme_style,
    "minimal" = theme_minimal,
    "classic" = theme_classic,
    "bw" = theme_bw,
    "light" = theme_light,
    theme_minimal
  )

  # Base plot
  p <- ggplot() +
    geom_point_rast(data = unSig, aes(x = lfc, y = logp), size = dotSize, colour = unSig.color) +
    geom_point_rast(data = u, aes(x = lfc, y = logp), col = up.color, alpha = alpha, size = dotSize) +
    geom_point_rast(data = d, aes(x = lfc, y = logp), col = down.color, alpha = alpha, size = dotSize) +
    labs(
      x = expression(log[2] ~ "Fold Change"),
      y = expression(-log[10] ~ "Adjusted P-value")
    ) +
    scale_x_continuous(limits = x_limits, breaks = lfc_breaks, expand = c(0, 0)) +
    theme_func(base_size = 12)

  # Y-axis scaling
  if (use_asinh) {
    p <- p + scale_y_continuous(trans = "asinh", limits = y_limits, breaks = padj_breaks, expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(limits = y_limits, breaks = padj_breaks, expand = c(0, 0))
  }

  # Count annotations
  if (count_position == "bottom") {
    p <- p + annotate("text", x = c(-Inf, Inf), y = c(-Inf, -Inf),
      hjust = c(-0.2, 1.2), vjust = c(-1.5, -1.5),
      label = c(paste("Down:", nrow(d)), paste("Up:", nrow(u))), size = 4, fontface = "bold")
  } else if (count_position == "top") {
    p <- p + annotate("text", x = c(-Inf, Inf), y = c(Inf, Inf),
      hjust = c(-0.2, 1.2), vjust = c(1.5, 1.5),
      label = c(paste("Down:", nrow(d)), paste("Up:", nrow(u))), size = 4, fontface = "bold")
  } else if (count_position == "corner") {
    p <- p + annotate("text", x = c(-Inf, Inf), y = c(Inf, Inf),
      hjust = c(-0.5, 1.5), vjust = c(2, 2),
      label = c(nrow(d), nrow(u)), size = 6)
  }

  # Load TF table if needed
  if (label_TF_only) {
    tf_table <- read.table(tf_table_path, header = TRUE, sep = "\t", quote = "")
  }

  # Select genes to label
  if (label_TF_only) {
    top_up <- u %>%
      mutate(gene = rownames(.)) %>%
      filter(gene %in% tf_table$Symbol) %>%
      arrange(desc(lfc)) %>%
      head(n_labels_up)
    top_down <- d %>%
      mutate(gene = rownames(.)) %>%
      filter(gene %in% tf_table$Symbol) %>%
      arrange(lfc) %>%
      head(n_labels_down)
  } else {
    top_up <- u %>%
      mutate(gene = rownames(.)) %>%
      arrange(desc(lfc)) %>%
      head(n_labels_up)
    top_down <- d %>%
      mutate(gene = rownames(.)) %>%
      arrange(lfc) %>%
      head(n_labels_down)
  }

  top_labels <- rbind(top_up, top_down)

  # Add labeled points and text
  if (nrow(top_labels) > 0) {
    p <- p + geom_point(
      data = top_labels,
      aes(x = lfc, y = logp),
      size = dotSize * 2,
      shape = 21,
      fill = ifelse(top_labels$lfc > 0, up.color, down.color),
      color = "black",
      stroke = 0.5
    )

    if (align_labels) {
      if (nrow(top_up) > 0) {
        p <- p + geom_text_repel(
          data = top_up,
          aes(x = lfc, y = logp, label = gene),
          size = 3, max.overlaps = Inf,
          box.padding = 0.4,
          segment.color = "grey50", segment.size = 0.3,
          hjust = 0, direction = "y",
          xlim = c(max(df$lfc) * 1.05, NA)
        )
      }
      if (nrow(top_down) > 0) {
        p <- p + geom_text_repel(
          data = top_down,
          aes(x = lfc, y = logp, label = gene),
          size = 3, max.overlaps = Inf,
          box.padding = 0.4,
          segment.color = "grey50", segment.size = 0.3,
          hjust = 1, direction = "y",
          xlim = c(NA, min(df$lfc) * 1.05)
        )
      }
	  } else {
      p <- p + geom_text_repel(
        data = top_labels,
        aes(x = lfc, y = logp, label = gene),
        size = 3, max.overlaps = 20,
        box.padding = 0.6,
        point.padding = 0.5,
        segment.color = "grey50", segment.size = 0.3,
        force = 2
      )
    }
  }

  return(p)
}