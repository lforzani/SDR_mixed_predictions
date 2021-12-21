rm(list=ls())

library(readr)
library(fastDummies)
library(MLmetrics)
library(pROC)
library(PCAmixdata)
library(MASS)

source("~/internal/PCAmix.cv.R")

#########LEER LOS DATOS

setone = read.csv("~/Data-Analysis/data-Krzanowski/setone.csv", header=FALSE, sep=";")

settwo = read.csv("~/Data-Analysis/data-Krzanowski/settwo.csv", header=FALSE, sep=";")

setthree = read.csv("~/Data-Analysis/data-Krzanowski/setthree.csv", header=FALSE, sep=";")

setfour = read.csv("~/Data-Analysis/data-Krzanowski/setfour.csv", header=FALSE, sep=";")


###################
# set one
setone = cbind(setone[,c(1,2,3,4,5,7)], log(setone[,6]), setone[,8:11])
datos = setone[,c(11,1:10)]
o=NULL
c=c(2:8)
b=c(9:11)

# set two
settwo = cbind(settwo[,1:3], log(settwo[,4]), settwo[,5:8])
datos=settwo[,c(8,1:7)]
o=NULL
c=c(2:5)
b=c(6:8)

# set three
setthree = cbind(setthree[,c(2,7,8)], log(setthree[,c(1,3,4,5,6)]), setthree[,9:13])
datos = setthree[,c(13,1:12)]
o=NULL
c=c(2:8)
b=c(9:13)

# set four
setfour = cbind(setfour[,1:4], log(setfour[,5:6]), setfour[,7:12])
datos= setfour[,c(12,1:11)]
o=NULL
c=c(2:7)
b=c(8:10)

################################

attach(datos)

n = dim(datos)[1]
K = n

aux.PCAmix = PCAmix.cv(datos, c, b, K)
MSE.PCAmix = aux.PCAmix$MSE
AUC.PCAmix = aux.PCAmix$AUC

aux.PCA = PCAmix.cv(datos, c, b = NULL, K)
MSE.PCA = aux.PCA$MSE
AUC.PCA = aux.PCA$AUC

round(c(MSE.PCA,MSE.PCAmix), 3)
round(c(AUC.PCA, AUC.PCAmix), 3)
