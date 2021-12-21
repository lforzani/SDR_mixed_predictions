#The following script is used to obatin results for Table 5 of the manuscript
#Sufficient Reductions in Regression with Mixed Predictors" Bura et al. (2021)

library(readr)
library(fastDummies)
library(MLmetrics)
library(PCAmixdata)
library(npreg)
library(np)

set.seed(1234)
#Read Data
aux <- read.csv("GovernanceDataExtended2.csv", header=TRUE, sep=";")

datos = aux[,c(1:7,10:20)]
attach(datos)
p = 6
q=11
o=NULL
Y=datos[1]
X=datos[2:7]
H=datos[8:18]
N0=100 

n = dim(datos)[1]


#read proyections (CG indices) computed using Matlab
PjPFC <- read.csv("ProyPFC.csv", header=FALSE, sep=",")
PjOptimal <- read.csv("ProyOptimal.csv", header=FALSE, sep=",")
PjSuboptimal1 <- read.csv("ProySuboptimal1.csv", header=FALSE, sep=",")
PjSuboptimal2 <- read.csv("ProySuboptimal2.csv", header=FALSE, sep=",")


##dimmension of the predictores
#p_dim = length(b)+length(c)


costos_PCA_MIX=matrix(0,nrow=n,ncol=1)

costos_PCA=matrix(0,nrow=n,ncol=1)

costos_NPPCA_MIX=matrix(0,nrow=n,ncol=1)

costos_NPPCA=matrix(0,nrow=n,ncol=1)

costos_PCA_factor=matrix(0,nrow=n,ncol=1)

costos_NPPFC=matrix(0,nrow=n,ncol=1)

costos_NPOptimal=matrix(0,nrow=n,ncol=1)

costos_NPSuboptimal=matrix(0,nrow=n,ncol=1)

for (i in 1:n) {
  
  Y_testing<-Y[i,]
  X_testing <-X[i,]
  H_testing <-H[i,]
 
  Y_training<-Y[-i,]
  X_training <-X[-i,]
  H_training <-H[-i,]
 
  H_testing_f<-apply(cbind(H_testing),2,factor)
  H_training_f<-apply(cbind(H_training),2,factor)
  
   ##############In M$score y M_PCA$score are linear combinations
  
  ## PCA mix
  M=PCAmix(X.quanti =  X_training, X.quali =H_training_f, ndim = 2, rename.level = TRUE,
           weight.col.quanti = NULL, weight.col.quali = NULL, graph = FALSE)
  
  ##Standard PCA
  M_PCA=prcomp( X_training,rank=1)
  
  
  M_scores=M$scores[,1]

  M_PCA_scores=predict(M_PCA,X_training)[,1]
  
  #####PCA_mix
  
  datos_PCA_mix=data.frame(XX=M_scores,Y=Y_training)
  
  lmfit_PCA_mix <- lm(Y~XX, data=datos_PCA_mix)
  
  bwnp_PCA_mix<-npregbw(Y~XX, data=datos_PCA_mix)
  
  np_PCA_mix <- npreg(bwnp_PCA_mix)
  
  #####PCA_standard
  datos_PCA=data.frame(XX=M_PCA_scores,Y=Y_training)
  
  lmfit_PCA <- lm(Y~XX, data=datos_PCA)
  
  bwnp_PCA<-npregbw(Y~XX, data=datos_PCA)
  
  np_PCA <- npreg(bwnp_PCA)
  
  #####PCA_standard with factors
  all=cbind(M_PCA_scores,H_training[,2:9])
  
  datos_factor= data.frame(XX=all, Y=Y_training)
  
  lmfit_PCA_factor <- lm(Y~XX.M_PCA_scores+ XX.BOL +XX.BRA +XX.CHL +XX.COL +XX.ECU +XX.PER +XX.PRY +XX.URY , data=datos_factor)
  
  
  
  ####Superivised reductions (reduction computed in Matlab)
  ##PFC
  pjPFCtr= PjPFC[2:n,i]
 
  datos_PFC=data.frame(XX=pjPFCtr,Y=Y_training)
  
  bwnp_PFC<-npregbw(Y~XX, data=datos_PFC)
 
  np_PFC <- npreg(bwnp_PFC)
  
  ###Optimal SDR for mixed predictors
 
   pjOptimaltr= PjOptimal[2:n,i]
  
   datos_Optimal=data.frame(XX=pjOptimaltr,Y=Y_training)
  
   bwnp_Optimal<-npregbw(Y~XX, data=datos_Optimal)
  
   np_Optimal <- npreg(bwnp_Optimal)
  
  ###Suboptimal SDR for mixed predictors
 
  datos_Suboptimal=data.frame(X1=PjSuboptimal1[2:n,i], X2=PjSuboptimal2[2:n,i],Y=Y_training)
  
  bwnp_Suboptimal<-npregbw(formula=Y~X1 + X2, data=datos_Suboptimal)
  
  np_Suboptimal <- npreg(bwnp_Suboptimal)
  
  #############################
  X_final_test=X_testing
  Cual_final_test=H_testing_f
  aux = names(Cual_final_test)
  
  #if the dimension of the conitnuos is 1 i need a matrix instead of a vector
  if(!is.null(X_final_test) & is.null(dim(X_final_test)[1])){X_final_test=as.matrix(X_final_test)}
  if(!is.null(Cual_final_test)){Cual_final_test=matrix(Cual_final_test, nrow = 1)
  colnames(Cual_final_test) <- aux}
  
  #predict on MIX (new scores)
  
  Mnew=predict(M,X.quanti=X_final_test,X.quali=Cual_final_test)
  
  

  #predict of PCA
  
  Mnew_PCA=predict(M_PCA,X_testing)
  
  
  ########################## 
  
  ##########B) prediction on the new score
  
  # PCAmix
  Mnew_datos = data.frame(XX=Mnew[,1])
  
  yhat_PCA_mix = predict(lmfit_PCA_mix, newdata =  Mnew_datos, type = "response") 
  nphat_PCA_mix = predict(np_PCA_mix, newdata=Mnew_datos, type="response")
  
  costos_PCA_MIX[i] = mean((yhat_PCA_mix - Y_testing)^2)
  costos_NPPCA_MIX[i] = mean((nphat_PCA_mix - Y_testing)^2)
  
  # PCA
  Mnew_datos = data.frame(XX=Mnew_PCA[,1])
  
  yhat_PCA = predict(lmfit_PCA, newdata =  Mnew_datos, type = "response") 
  nphat_PCA = predict(np_PCA, newdata=Mnew_datos, type="response")
  
  costos_PCA[i] = mean((yhat_PCA - Y_testing)^2)
  costos_NPPCA[i] = mean((nphat_PCA - Y_testing)^2)
  
# PCA Factor
  
  M_PCA_scores=Mnew_PCA[,1]
  auxtest=cbind(M_PCA_scores,H_testing[,2:9])
  Mnew_datos_factor = data.frame(XX=auxtest)
  yhat_PCA_factor = predict(lmfit_PCA_factor, newdata = Mnew_datos_factor, type = "response") 
  
  costos_PCA_factor[i] = mean((yhat_PCA_factor - Y_testing)^2)
  
  #PFC
  Mnew_datos = data.frame(XX=PjPFC[1,i])
  nphat_PFC = predict(np_PFC, newdata=Mnew_datos, type="response")
  
  costos_NPPFC[i] = mean((nphat_PFC - Y_testing)^2)
  
  #SDR mix optimal
  Mnew_datos = data.frame(XX=PjOptimal[1,i])
  
  nphat_Optimal = predict(np_Optimal, newdata=Mnew_datos, type="response")
  
  costos_NPOptimal[i] = mean((nphat_Optimal - Y_testing)^2)
  
  #SDR mix Suboptimal
  Mnew_datos = data.frame(X1=PjSuboptimal1[1,i], X2=PjSuboptimal2[1,i])
  
  nphat_Suboptimal = predict(np_Suboptimal, newdata=Mnew_datos, type="response")
  
  costos_NPSuboptimal[i] = mean((nphat_Suboptimal - Y_testing)^2)
  
  
  }

