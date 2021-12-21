rm(list=ls())
library(readr)

n = c(100, 200, 300, 500, 750)

########################################## FIGURE 2 ################################################

# Dim 1
Continuous_d1_estimator <- read_csv("/Figures/MixedPredictors/Mixed_estimator_dim1_suboptimal_cont.csv", col_names = FALSE)

Continuous_d1_projection <- read_csv("/Figures/MixedPredictors/Mixed_proj_dim1_suboptimal_cont.csv", col_names = FALSE)

Binaria_d1_estimator <- read_csv("/Figures/MixedPredictors/Mixed_estimator_dim1_suboptimal_disc.csv", col_names = FALSE)

Binaria_d1_projection <- read_csv("/Figures/MixedPredictors/Mixed_proj_dim1_suboptimal_disc.csv", col_names = FALSE)


# Dim 2
Continuous_d2_estimator <- read_csv("/Figures/MixedPredictors/Mixed_estimator_dim2_suboptimal_cont.csv", col_names = FALSE)

Continuous_d2_projection <- read_csv("/Figures/MixedPredictors/Mixed_estimator_dim2_suboptimal_disc.csv", col_names = FALSE)

Binaria_d2_estimator <- read_csv("/Figures/MixedPredictors/Mixed_proj_dim2_suboptimal_cont.csv", col_names = FALSE)

Binaria_d2_projection <- read_csv("/Figures/MixedPredictors/Mixed_proj_dim2_suboptimal_disc.csv", col_names = FALSE)

################

dfpd1c <- data.frame(meanX1  = colMeans(Continuous_d1_estimator, na.rm = TRUE), meanPX1 = colMeans(Continuous_d1_projection, na.rm = TRUE))
dfpd1b <- data.frame(meanX1  = colMeans(Binaria_d1_estimator, na.rm = TRUE), meanPX1 = colMeans(Binaria_d1_projection, na.rm = TRUE))

dfpd2c <- data.frame(meanX1  = colMeans(Continuous_d2_estimator, na.rm = TRUE), meanPX1 = colMeans(Continuous_d2_projection, na.rm = TRUE))
dfpd2b <- data.frame(meanX1  = colMeans(Binaria_d2_estimator, na.rm = TRUE), meanPX1 = colMeans(Binaria_d2_projection, na.rm = TRUE))


# Start Plot --------------------

# 1. Open eps file
setEPS()
postscript("Figure-2.eps")

# 2. Create the plot
m <- matrix(c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.5,0.5))

# Plots --------------------
plot(dfpd1c[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd1c[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n")
lines(dfpd1b[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd1b[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)

plot(dfpd2c[,1], type = 'l', lty = 1, lwd=1.5, cex.axis = 1, cex.main = 1, main = '', ylim = c(0,1), cex.lab = 1, ylab = "Error", xlab = 'Sample Size', xaxt = 'n')
axis(side = 1,  tick = TRUE, at = 1:length(n), labels = n, cex.axis = 1, cex.lab = 1)
lines(dfpd2c[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n") ##############
lines(dfpd2b[,1], lty = 1, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)
lines(dfpd2b[,2], lty = 2, lwd=1.5, main = '', ylim = c(0,1), xaxt="n", col = 2)

# legend ----------------------
legend.name = c('Continuous Suboptimal Reduction', 'Continuous Suboptimal Projection', 'Binary Suboptimal Reduction', 'Binary Suboptimal Projection') 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",  bty = "n", xpd = TRUE, horiz = FALSE, xjust = 0, yjust = 0, legend = legend.name, cex = 1, lty = 1:2, lwd=1.5, title = "", col= c(1,1,2,2,3,3))


# 3. Close the file
dev.off()
