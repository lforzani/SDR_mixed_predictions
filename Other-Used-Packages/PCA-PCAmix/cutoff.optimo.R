cutoff.optimo = function(Y, phat){
  # given the truth Y and the fitted probabilities phat
  # this function computes the optimal cutoff minimizing the IU 
  
  library(pROC)
  
  aa = roc(Y, phat,quiet=TRUE)

  cutoffs = aa$thresholds
  AUC = aa$auc
  TPR = aa$sensitivities
  TNR = aa$specificities
  IU = abs(TPR-AUC) + abs(TNR-AUC)
  
  #plot(IU, type = 'l') 
  
  cutoff.optimo = cutoffs[which.min(IU)]
  valor = IU[which.min(IU)]
  
  return(c(cutoff.optimo, valor))
}

