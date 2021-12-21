rm(list=ls())
library(readr)

n = c(100, 200, 300, 500, 750)

########################################## FIGURE 1 ################################################

# Start Plot

# 1. Open eps file
setEPS()
postscript("Figure-1.eps")

# 2. Create the plot
m <- matrix(c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = m)

#################################### Continuous Predictors #########################################

# Read Data ---------------
# Dim 1
Continuous_d1_estimador <- read_csv("~/Figures/ContinuousPredictors/Continuous-d1-estimator.csv", col_names = FALSE)

Continuous_d1_projection <- read_csv("~/Figures/ContinuousPredictors/Continuous-d1-projection.csv", col_names = FALSE)

Continuous_d1_estimadorPCA <- read_csv("~/Figures/ContinuousPredictors/Continuous-d1-estimator_PCA.csv", col_names = FALSE)

Continuous_d1_projection_PCA <- read_csv("~/Figures/ContinuousPredictors/Continuous-d1-projection_PCA.csv", col_names = FALSE)


# Dim 2
Continuous_d2_estimador <- read_csv("~/Figures/ContinuousPredictors/Continuous-d2-estimator.csv", col_names = FALSE)

Continuous_d2_projection <- read_csv("~/Figures/ContinuousPredictors/Continuous-d2-projection.csv", col_names = FALSE)

Continuous_d2_estimadorPCA <- read_csv("~/Figures/ContinuousPredictors/Continuous-d2-estimator_PCA.csv", col_names = FALSE)

Continuous_d2_projection_PCA <- read_csv("~/Figures/ContinuousPredictors/Continuous-d2-projection_PCA.csv", col_names = FALSE)


# Means --------------------

# Optimal
dfpd1 <- data.frame(meanX1  = colMeans(Continuous_d1_estimador), meanPX1  = colMeans(Continuous_d1_projection))
dfpd2 <- data.frame(meanX1  = colMeans(Continuous_d2_estimador), meanPX1  = colMeans(Continuous_d2_projection ))

# PCA
dfpd1PCA <- data.frame(meanX1  = colMeans(Continuous_d1_estimadorPCA), meanPX1  = colMeans(Continuous_d1_projection_PCA))
dfpd2PCA <- data.frame(meanX1  = colMeans(Continuous_d2_estimadorPCA), meanPX1  = colMeans(Continuous_d2_projection_PCA ))

### Plots 
plot(dfpd1[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = 'd = 1', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd1[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n")
lines(dfpd1PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)

plot(dfpd2[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = 'd = 2', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd2[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1))
lines(dfpd2PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)


#################################### Binary Predictors #########################################

# Read Data 

# Dim 1
Binaria_d1_estimador <- read_csv("~/Figures/BinaryPredictors/Binary-d1-estimator.csv", col_names = FALSE)

Binaria_d1_projection <- read_csv("~/Figures/BinaryPredictors/Binary-d1-projection.csv", col_names = FALSE)

Binaria_d1_estimadorPCA <- read_csv("~/Figures/BinaryPredictors/Binary-d1-estimator_PCA.csv", col_names = FALSE)

Binaria_d1_projection_PCA <- read_csv("~/Figures/BinaryPredictors/Binary-d1-estimator_PCA.csv", col_names = FALSE)

Binaria_d1_estimador_PFC <- read_csv("~/Figures/BinaryPredictors/Binary-d1-estimator_PFC.csv", col_names = FALSE)

Binaria_d1_estimador_PFC <- read_csv("~/Figures/BinaryPredictors/Binary-d1-estimator_PFC.csv", col_names = FALSE)

# Dim 2
Binaria_d2_estimador <- read_csv("~/Figures/BinaryPredictors/Binary-d2-estimator.csv", col_names = FALSE)

Binaria_d2_projection <- read_csv("~/Figures/BinaryPredictors/Binary-d2-projection.csv", col_names = FALSE)

Binaria_d2_estimadorPCA <- read_csv("~/Figures/BinaryPredictors/Binary-d2-estimator_PCA.csv", col_names = FALSE)

Binaria_d2_projection_PCA <- read_csv("~/Figures/BinaryPredictors/Binary-d2-projection_PCA.csv", col_names = FALSE)

Binaria_d2_estimador_PFC <- read_csv("~/Figures/BinaryPredictors/Binary-d2-estimator_PFC.csv", col_names = FALSE)

Binaria_d2_estimador_PFC <- read_csv("~/Figures/BinaryPredictors/Binary-d2-projection_PFC.csv", col_names = FALSE)


# Means --------------------

# Optimal
dfpd1 <- data.frame(meanX1  = colMeans(Binaria_d1_estimador, na.rm = TRUE), meanPX1 = colMeans(Binaria_d1_projection, na.rm = TRUE))
dfpd2 <- data.frame(meanX1  = colMeans(Binaria_d2_estimador, na.rm = TRUE), meanPX1 = colMeans(Binaria_d2_projection, na.rm = TRUE))

# PCA
dfpd1PCA <- data.frame(meanX1  = colMeans(Binaria_d1_estimadorPCA, na.rm = TRUE), meanPX1  = colMeans(Binaria_d1_projection_PCA, na.rm = TRUE))
dfpd2PCA <- data.frame(meanX1  = colMeans(Binaria_d2_estimadorPCA, na.rm = TRUE), meanPX1  = colMeans(Binaria_d2_projection_PCA, na.rm = TRUE))

# PFC
dfpd1PFC <- data.frame(meanX1  = colMeans(Binaria_d1_estimador_PFC, na.rm = TRUE), meanPX1  = colMeans(Binaria_d1_estimador_PFC, na.rm = TRUE))
dfpd2PFC <- data.frame(meanX1  = colMeans(Binaria_d2_estimador_PFC, na.rm = TRUE), meanPX1  = colMeans(Binaria_d2_estimador_PFC, na.rm = TRUE))

## Plots
plot(dfpd1[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd1[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n")
lines(dfpd1PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1PFC[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)
lines(dfpd1PFC[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)

plot(dfpd2[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd2[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1))
lines(dfpd2PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2PFC[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)
lines(dfpd2PFC[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)
#################################### Mixed Predictors #########################################

# Read Data

# Dim 1
Mixed_d1_estimador <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim1_optimal.csv", col_names = FALSE)

Mixed_d1_projection <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim1_optimal.csv", col_names = FALSE)

Mixed_d1_estimadorPCA <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim1_PCA.csv", col_names = FALSE)

Mixed_d1_projection_PCA <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim1_PCA.csv", col_names = FALSE)

Mixed_d1_estimador_PFC <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim1_PFC.csv",col_names = FALSE)

Mixed_d1_estimador_PFC <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim1_PFC.csv", col_names = FALSE)


# Dim 2
Mixed_d2_estimador <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim2_optimal.csv", col_names = FALSE)

Mixed_d2_projection <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim2_optimal.csv", col_names = FALSE)

Mixed_d2_estimadorPCA <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim2_PCA.csv", col_names = FALSE)

Mixed_d2_projection_PCA <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim2_PCA.csv", col_names = FALSE)

Mixed_d2_estimador_PFC <- read_csv("~/Figures/MixedPredictors/Mixed_estimator_dim2_PFC.csv", col_names = FALSE)

Mixed_d2_estimador_PFC <- read_csv("~/Figures/MixedPredictors/Mixed_proy_dim2_PFC.csv", col_names = FALSE)


# Means --------------------

# Optimal
dfpd1 <- data.frame(meanX1  = colMeans(Mixed_d1_estimador, na.rm = TRUE), meanPX1 = colMeans(Mixed_d1_projection, na.rm = TRUE))
dfpd2 <- data.frame(meanX1  = colMeans(Mixed_d2_estimador, na.rm = TRUE), meanPX1 = colMeans(Mixed_d2_projection, na.rm = TRUE))

# PCA
dfpd1PCA <- data.frame(meanX1  = colMeans(Mixed_d1_estimadorPCA, na.rm = TRUE), meanPX1  = colMeans(Mixed_d1_projection_PCA, na.rm = TRUE))
dfpd2PCA <- data.frame(meanX1  = colMeans(Mixed_d2_estimadorPCA, na.rm = TRUE), meanPX1  = colMeans(Mixed_d2_projection_PCA, na.rm = TRUE))

# PFC
dfpd1PFC <- data.frame(meanX1  = colMeans(Mixed_d1_estimador_PFC, na.rm = TRUE), meanPX1  = colMeans(Mixed_d1_estimador_PFC, na.rm = TRUE))
dfpd2PFC <- data.frame(meanX1  = colMeans(Mixed_d2_estimador_PFC, na.rm = TRUE), meanPX1  = colMeans(Mixed_d2_estimador_PFC, na.rm = TRUE))

## Plots
plot(dfpd1[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd1[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n")
lines(dfpd1PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1PFC[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)
lines(dfpd1PFC[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)

plot(dfpd2[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd2[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1))
lines(dfpd2PCA[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2PCA[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2PFC[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)
lines(dfpd2PFC[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 3)


## Legend 
legend.name = c('Optimal Reduction', 'Optimal Projection','PCA Reduction', 'PCA Projection', 'PFC Reduction', 'PFC Projection') 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", bty = "n", xpd = TRUE, horiz = FALSE, xjust = 0, yjust = 0, legend = legend.name, cex = 1, lty = 1:2, lwd=1.5, title = "", col= c(1,1,2,2,3,3))

# 3. Close the file
dev.off()

