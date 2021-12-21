PCAmix.cv = function(datos, c, b, K){
  # if K = k performs k-fold cv
  # if K = n performs leave-one-out cv
  
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/Datos/data-Krzanowski/cutoff.optimo.R")
  
  library(readr)
  library(fastDummies)
  library(MLmetrics)
  library(pROC)
  library(PCAmixdata)
  library(MASS)
  
  n = dim(datos)[1]
  
  error = NULL
  auc = NULL
  Yhat = NULL
  
  folds <- cut(sample(seq(1, nrow(datos))), breaks=K, labels=FALSE)
  
  for (i in 1:K){
    cat(" i = ", i)
    
    if (K == n){
      idx <-  i
    } else {
      idx <- which(folds == i, arr.ind=TRUE)
    }

    testing <- datos[idx,]
    training <- datos[-idx,]
    
    # train -------------------------
    if (is.null(b)){
      X.quali = NULL
      X.quanti = training[,c]
    } else {
      X.quali = apply(cbind(training[,b]), 2, factor)
      X.quanti = training[,c]
    }
   
    PCAmix.obj = PCAmix(X.quanti = X.quanti, X.quali = X.quali, ndim = 2, rename.level = TRUE, graph = FALSE)
    
    # saco los coeficientes de la training (para usarlos despues) 
    #alfas = PCAmix.obj$coef$dim1[,1]
    
    # scores
    M.train = PCAmix.obj$scores[,1]
    
    # nuevos datos
    data.new <- data.frame(X = M.train, Y = training[,1])
    
    # ajusto un glm de y~scores
    fit <- glm(Y ~ X, data = data.new, family = 'binomial'(link = "logit"))
      
    corte <- cutoff.optimo(training[,1], fit$fitted.values)[1] 
    
    # test -------------------------
    if (is.null(b)){
      X.quali = NULL
      X.quanti = testing[,c]
      #X = cbind(X.quanti)
    } else {
      # usando predict
      X.quanti = testing[,c]
      X.quali = apply(cbind(testing[,b]), 2, factor) 
      aux = names(testing[,b]) #names(X.quali)
      # para que ande el predict
      if(!is.null(X.quanti) & is.null(dim(X.quanti)[1])){X.quanti=as.matrix(X.quanti)}
      if(!is.null(X.quali)){X.quali=matrix(X.quali, ncol = length(b), nrow = dim(testing)[1])
      colnames(X.quali) <- aux}
    }
    
    # scores en la testing
    M.test = predict(PCAmix.obj, X.quanti = X.quanti, X.quali = X.quali)[,1]
  
    # ajusto glm en la testing
    data.new <- data.frame(X = M.test)
    phat = predict(fit, newdata = data.new, type = "response")
  
    if (K == n){
      Yhat[i] <- 1*(phat > corte)
      error[i] = 1-mean(testing[,1]==Yhat[i])
    } else {
      Yhat <- 1*(phat > corte)
      error[i] = 1-mean(testing[,1]==Yhat)
      auc[i] =  roc(testing[,1], Yhat, quiet = TRUE)$auc}
  }
  
  MSE = mean(error)
  
  if (K == n){
    AUC = roc(datos[,1], Yhat, quiet = TRUE)$auc
  } else {
    AUC = mean(auc)
  }
  
  return(list(MSE = MSE, AUC = AUC, cutoff.optimo = corte))
}
