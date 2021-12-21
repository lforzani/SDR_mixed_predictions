#This script compute Governance Index using PCAmix method

library(MASS)
library(car)
library(PCAmixdata)

data <- read.csv("GovernanceDataExtended2.csv", header=TRUE, sep=";")
datos = data[,c(1:7,10:20)]
attach(datos)
p = 6
q=11
Y=datos[1]
X=datos[2:7]
H=datos[8:18]

n = dim(datos)[1]

Haux<-as.matrix(H)
Hf<-apply(cbind(Haux),2,factor)
pcamixto<-PCAmix(X.quanti = X, X.quali = Hf, ndim = 2,rename.level = TRUE, weight.col.quanti = NULL, weight.col.quali = NULL, graph = TRUE)
projPCAmix <-pcamixto$scores[,1]

CGindex_PCAmix = (projPCAmix-min(projPCAmix))/(max(projPCAmix)-min(projPCAmix))

#save PCA_mix
write.csv(CGindex_PCAmix, file = "PCAmixIndex.csv")
