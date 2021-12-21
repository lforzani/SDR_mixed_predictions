#library(MASS)
###inputs

#datos: is a matrix with columns conitnuous, binary and ordinal variables, 
#we give outside the index as c=indexes for conitnuos, o=indexes for ordinals
#and b=indexes for binary


#N0: since we are asking all the posible choise of pstar variables on p
#variable it can be too much and we limite to N0 the number of choises

#gamma: if 1 choose among all combination, smaller choose the best among some
#combinations.

#output: 
#we are looking for  pstar variables among the p that the PCA linear combinations
#are the closest to the PCA with all the variables

#ultimas_PCAmix: using PCAmix done by the pacakge PCAmixdata (reference)

#ultimas_PCA: using regular PCA where the ordinals are converted into dicotomics


PCA_ambos<-function(datos,N0,gamma,p,pstar){
  
  library(PCAmixdata)
  library(fastDummies)
  
  ##look for all the combination of pstar variables onto p. If the number 
  ##of combinations is too large it choose only N0 of them (randomly)
  
  cc = choose(p, pstar)
  if (cc < N0){
    MM = combn(2:(p+1), pstar)
    gamma = 1
  } else {
    MM = matrix(0, nrow = pstar, ncol = N0)
    for (i in 1:N0){
      MM[,i] = sample(2:(p+1),pstar,replace=FALSE) 
    }
  }
  nn = dim(MM)[2]
  
  result_PCA=NULL
  result_PCAmix=NULL
  subset=matrix(0,nrow=nn,ncol=pstar)
  
  ###PCA MIX
  #take the continuos variables
  X=cbind(datos[,c])
  #take the binary data
  H=datos[,b]
  
  
  # to be able to use PCAmix it needs to add together the binary and ordinal data
  if(is.null(o)){
    Cual=apply(cbind(datos[,b]), 2, factor)}
  if(!is.null(o)){
    Cual=apply(cbind(datos[,b],datos[,o]), 2, factor)}
  
  
  
  #apply PCAmix
  M=PCAmix(X.quanti = X, X.quali =Cual, ndim = 2, rename.level = TRUE,
           weight.col.quanti = NULL, weight.col.quali = NULL, graph = FALSE)
  
  ##get data*first PCA
  Morig_PCAmix=M$scores[,1]
  
  
  
  ##PCA 
  #convert the dummies into categoricals
  if(is.null(o)){
    XX=cbind(X,H)}
  if(!is.null(o)){
    O_dummy=dummy_cols(datos[,o],select_columns=c(names(datos[,o])))[,-(1:length(o))]
    XX=cbind(X,H,O_dummy)}
  
  #apply regular PCA
  M_comun=prcomp(XX)
  
  ##get data*first PCA
  Morig_PCA=predict(M_comun,XX)[,1]
  
  
  for (s in 1:nn){
    subset[s,] = MM[,s]
    
    
    ###accomodate data for applying PCA and PCAmix
    
    
    acomodate_data=accomodate_variables_for_PCA_PCAmix(subset[s,],c,b,o,datos)
    
    ### apply PCAmix to subset
    
    Xsubset=acomodate_data$Xsubset
    Cualsubset=acomodate_data$Cualsubset
    M=PCAmix(X.quanti = Xsubset, X.quali =Cualsubset, ndim = 2, rename.level = TRUE,
             weight.col.quanti = NULL, weight.col.quali = NULL, graph = FALSE)
    
    Mnew_PCAmix=M$scores[,1]
    
    ### apply PCA to subset   
    
    XXsubset = acomodate_data$XXsubset
    M_PCA=prcomp(XXsubset)
    
    
    Mnew_PCA=predict(M_PCA,XXsubset)[,1]
    
    
    #compute the diference between PCA and PCAmix
    
    result_PCAmix[s]=(t(Morig_PCAmix)%*%Morig_PCAmix)^(-1)*(t(Morig_PCAmix)%*%Mnew_PCAmix)%*%(t(Mnew_PCAmix)%*%Mnew_PCAmix)^(-1)%*%t(Mnew_PCAmix)%*%Morig_PCAmix
    
    result_PCA[s]=(t(Morig_PCA)%*%Morig_PCA)^(-1)*(t(Morig_PCA)%*%Mnew_PCA)%*%(t(Mnew_PCA)%*%Mnew_PCA)^(-1)%*%t(Mnew_PCA)%*%Morig_PCA
    
    result_PCAmix[s]=abs(result_PCAmix[s])
    result_PCA[s]=abs(result_PCA[s])
    
    
  }
  
  #find the subset that gives the closest linear combination
  
  ##for PCAmix
  best <- order(result_PCAmix, decreasing=TRUE)[1:(gamma*nn)]
  best_subset <- subset[best, ]
  
  freqTable <- table(best_subset)
  ultimas_PCAmix=as.numeric(names(freqTable)[order(freqTable, decreasing = TRUE)])
  
  ##for PCA
  
  best <- order(result_PCA, decreasing=TRUE)[1:(gamma*nn)]
  best_subset <- subset[best, ]
  
  freqTable <- table(best_subset)
  ultimas_PCA=as.numeric(names(freqTable)[order(freqTable, decreasing = TRUE)])
  
  
  
  return(rbind(ultimas_PCAmix[1:pstar],ultimas_PCA[1:pstar]))
  
  
}


