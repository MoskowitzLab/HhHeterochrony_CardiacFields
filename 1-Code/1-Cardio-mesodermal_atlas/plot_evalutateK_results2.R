plot_evalutateK_results2 <- function(aicMat, k = NULL, aicDiff = 2) {
  if (is.null(k)) {
    k <- as.numeric(sub("k: ", "", colnames(aicMat)))
    }
  
  init_shape <- graphics::par()$mfrow
  graphics::par(mfrow = c(1, 4))
  # boxplots of AIC
  devs <- matrix(NA, nrow = nrow(aicMat), ncol = length(k))
  for (ii in seq_len(nrow(aicMat))) {
    devs[ii, ] <- aicMat[ii, ] - mean(aicMat[ii, ])
  }
  graphics::boxplot(devs, ylab = "Deviation from genewise average AIC",
                    xlab = "Number of knots", xaxt = "n")
  graphics::axis(1, at = seq_len(length(k)), labels = k)
  # scatterplot of average AIC
  plot(x = k, y = colMeans(aicMat, na.rm = TRUE), type = "b",
       ylab = "Average AIC", xlab = "Number of knots")
  # scatterplot of relative AIC
  plot(x = k, y = colMeans(aicMat / aicMat[, 1], na.rm = TRUE), type = "b",
       ylab = "Relative AIC", xlab = "Number of knots")
  # barplot of optimal AIC for genes with at least a difference of 2.
  aicRange <- apply(apply(aicMat, 1, range), 2, diff)
  varID <- which(aicRange > aicDiff)
  if (length(varID) > 0) {
    aicMatSub <- aicMat[varID, ]
    tab <- table(k[apply(aicMatSub, 1, which.min)])
    graphics::barplot(tab, xlab = "Number of knots",
                      ylab = "# Genes with optimal k")
  }
  graphics::par(mfrow = init_shape)
  return()
}
