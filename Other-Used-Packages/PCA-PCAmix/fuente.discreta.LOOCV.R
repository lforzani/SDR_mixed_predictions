fuente.discreta.LOOCV = function(Base_mixed, o = NULL, c, b){
  # el orden en Base_mixed debe ser: Y, Continuas, Binarias
  # o, c, b son los indices donde estan las variables Ordinales (aca es nulo),
  # Continuas y Binarias
  library(readr)
  library(fastDummies)
  library(MLmetrics)
  
  trh = seq(0.2,0.8,0.05)
  
  ############## funciones que se necesitan ##############
  # path pame
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/funciones PCA-PCAmix/cutoff.optimo.minMSE.R")
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/funciones PCA-PCAmix/PCA_ambos.R")
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/funciones PCA-PCAmix/accomodate_variables_for_PCA_PCAmix.R")
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/funciones PCA-PCAmix/accomodate_PCA_2.R")
  source("~/Dropbox/INVESTIGACION/diego-pame-rodrigo (1)/Mixto/Codigo_Mixtos/independiente2019/funciones PCA-PCAmix/accomodate_PCA_MIX_2.R")
  
  #opcion 2
  N0=1000
  
  n = dim(Base_mixed)[1]
  
  ############## definition of variables #############
  # o: ordinal
  # b: binaria
  # c: continua
  
  
  
  ##dimmension of the predictores
  p_dim = length(o)+length(b)+length(c)
  
  
  ####################################
  ####################################
  ####################################
  ### 
  ###################################
  ######################################
  ###################################
  ####set to 0
  costos_PCA_MIX=matrix(0,nrow=n,ncol=p_dim)
  costos_PCA=matrix(0,nrow=n,ncol=p_dim)
  
  # MSE_PCA_mix=matrix(0,nrow=n,ncol=p_dim)
  # MSE_PCA=matrix(0,nrow=n,ncol=p_dim)
  
  #creat the 5 folders
  set.seed(1000)
  
  # we chose how many variables give the minimun error
  # of prediction
  
  #phat_PCA_MIX = NULL
  #phat_PCA = NULL
  
  #cat("== starting ps procedure for PFCMix (this can take several minutes)===")
  for (ps in 1:p_dim){
    #cat("finished p0 = ", ps)
    print(ps)
    ######### elegimos el corte optimo #######
  
    
    for (i in 1:n) {
    
     # cat(", finished i = ", i)
      #Randomly shuffle the data
      ######A) datos repartii=1dos   
      #Segement your data by fold using the which() function 
      M_l_testing <- Base_mixed[i, ]
      M_l_training <- Base_mixed[-i, ]
      
      ############################
      ##############Working on training
      ############################
      ### choosing the ps indexes that gives the combinations closes to PCA and PCA mix
      #######using all the variables
      
      indices=PCA_ambos(M_l_training,N0,.1,p_dim,ps)
      
      
      ####### Find the linear combinations, for that we need to acoomodate the variables 
      
      
      #####acomodo los indices para PCA comun PCAmix
      
      acomodate_data=accomodate_variables_for_PCA_PCAmix(indices[2,],c,b,o,M_l_training)
      
      
      ###PCA comun
      Total_final_PCA=acomodate_data$XXsubset
      
      
      ###PCA mix
      
      acomodate_data=accomodate_variables_for_PCA_PCAmix(indices[1,],c,b,o,M_l_training)
      
      X_final=acomodate_data$Xsubset
      Cual_final=acomodate_data$Cualsubset
      
      ##############
      ##############en M$score y M_PCA$score estan las combinaciones lineales
      
      
      ##mix
      M=PCAmix(X.quanti = X_final, X.quali =Cual_final, ndim = 2, rename.level = TRUE,
               weight.col.quanti = NULL, weight.col.quali = NULL, graph = FALSE)
      
      
      
      ##pc comun
      M_PCA=prcomp(Total_final_PCA,rank=1)
      
      
      ##########C) aplico glmfit y glm a los scorse (combinaciones lineales de PCA) en el training        
      
      #####miro los scores de cada una en el training
      M_scores=M$scores[,1]
      M_PCA_scores=predict(M_PCA,Total_final_PCA)[,1]
      
      #####PCA_mix
      datos_PCA_mix=data.frame(XX=M_scores,Y=M_l_training[,1])
      glmfit_PCA_mix <- glm(Y~XX, data=datos_PCA_mix, family = 'binomial'(link = "logit"))
      
      cutoff_optimo_pca_mix = cutoff.optimo.minMSE(glmfit_PCA_mix,M_l_training[,1],M_scores,trh)
      
      #####PCA_comun
      datos_PCA=data.frame(XX=M_PCA_scores,Y=M_l_training[,1])
      glmfit_PCA <- glm(Y~XX, data=datos_PCA, family = 'binomial'(link = "logit"))
      
      cutoff_optimo_pca = cutoff.optimo.minMSE(glmfit_PCA,M_l_training[,1],M_PCA_scores,trh) 
      
      ##############
      ##############
      ##############Working on the test
      ##############
      ###############
      #######A) see the liPampeanar combination on the test
      
      #####acomodate the index for PCA and PCAmix
      
      acomodate_data=accomodate_variables_for_PCA_PCAmix(indices[2,],c,b,o,M_l_testing)
      
      ###PCA comun
      Total_final_PCA_test=acomodate_data$XXsubset
      
      ###PCA mix
      acomodate_data=accomodate_variables_for_PCA_PCAmix(indices[1,],c,b,o,M_l_testing)
      
      X_final_test=acomodate_data$Xsubset
      Cual_final_test=acomodate_data$Cualsubset
      aux = names(Cual_final_test)
      
      #if the dimension of the conitnuos is 1 i need a matrix instead of a vector
      if(!is.null(X_final_test) & is.null(dim(X_final_test)[1])){X_final_test=as.matrix(X_final_test)}
      if(!is.null(Cual_final_test)){Cual_final_test=matrix(Cual_final_test, nrow = 1)
      colnames(Cual_final_test) <- aux}
      
      
      #############################
      #predict on MIX (new scores)
      
      Mnew=predict(M,X.quanti=X_final_test,X.quali=Cual_final_test)
      
      
      aux=setdiff(names(Total_final_PCA),names(Total_final_PCA_test))
      if (length(aux)!=0){
        Total_final_PCA_test[aux]=0}
      
      #predict of PCA
      
      Mnew_PCA=predict(M_PCA,Total_final_PCA_test)
      
      
      ########################## 
      
      ##########B) predcition on the new score
      
      # PCAmix
      Mnew_datos = data.frame(XX=Mnew[,1],Y=M_l_testing[,1])
      
      Mnew_aux_PCA_mix_glm = predict(glmfit_PCA_mix, newdata =  Mnew_datos, type = "response") 
      
      yhat_PCA_mix = 1*(Mnew_aux_PCA_mix_glm > cutoff_optimo_pca_mix)
      
      costos_PCA_MIX[i,ps] = 1-mean(yhat_PCA_mix == M_l_testing[,1])
    
     
       # PCA
       Mnew_datos = data.frame(XX=Mnew_PCA[,1],Y=M_l_testing[,1])
      
       Mnew_aux_PCA_glm = predict(glmfit_PCA, newdata =  Mnew_datos, type = "response") 
       
       yhat_PCA = 1*(Mnew_aux_PCA_glm > cutoff_optimo_pca)
       
       costos_PCA[i,ps] = 1-mean(yhat_PCA == M_l_testing[,1])
      
       
       #########C) save the predictions 
      
       # resultados_PCA_MIX_glm = data.frame("fiteados"= Mnew_aux_PCA_mix_glm, "verdaderos"=M_l_testing[,1])
       # 
       # 
       # resultados_PCA_glm = data.frame("fiteados"= Mnew_aux_PCA_glm, "verdaderos"=M_l_testing[,1])
       # 
      
      #######D) miro para cada test i la funcion que estoy mirando que hace lo mejor posible
      
      #phat_PCA_MIX[i] = resultados_PCA_MIX_glm$fiteados
      #phat_PCA[i] = resultados_PCA_MIX_glm$fiteados
  
      
      #Yverdaderas = Base_mixed[1,]
      
      # guardar los modelos para usarlos ocon el cutoff optimo
      # para clasificar 
      
      # aux_pca_mix = cutoff.optimo(Yverdaderas, phat_PCA_MIX) 
      # costos_PCA_MIX[i,ps]=aux_pca_mix[2]
      # aux_pca = cutoff.optimo(Yverdaderas, phat_PCA) 
      # costos_PCA[i,ps]=aux_pca[2]
    
    }
  }
  
  #looking for the ps in each cases that makes minimum cost 
  ps_optimo_PCA_glm=which.min(apply(costos_PCA,2,mean))
  
  ps_optimo_PCA_mix_glm=which.min(apply(costos_PCA_MIX,2,mean))
  
  
  #######Now, I have to do everything again since I did not save the variables for those ps. I have to
  ######to crossvalidation but for the same folders.
  
  
 # cat("== starting CV procedure for PFCMix (this can take several minutes)===")
  
  ###########Second part: finding the variables.
  ############
  ############
  ############
  ############
  ############
  ############
  ############
  ############
  
  #in each folder we found the best p_optimas variables, laeter we see which one appear more.
  
  indice_final_PCA_glm=matrix(0,nrow=n,ncol=ps_optimo_PCA_glm)
  indice_final_PCA_MIX_glm=matrix(0,nrow=n,ncol=ps_optimo_PCA_mix_glm)
  
  # MSE_PCA=NULL
  # MSE_PCA_mix=NULL
  # 
  # MSE_PCA_T=NULL
  # MSE_PCA_mix_T=NULL
  
  costos_PCA=NULL
  costos_PCA_MIX=NULL
 
  # corte_PCA_mix=NULL
  # corte_PCA=NULL
  
  # corte_PCA_mix_T=NULL
  # corte_PCA_T=NULL
  # 
  # 
  # tabla_pca_mix_glm=matrix(NA,ncol=4,nrow=5)
  # tabla_pca_glm=matrix(NA,ncol=4,nrow=5)
  # 
  # tabla_pca_mix_glm_T=matrix(NA,ncol=4,nrow=5)
  # tabla_pca_glm_T=matrix(NA,ncol=4,nrow=5)
  
  for (i in 1:n){
    print(i)
    #cat(", finished i = ", i)

    ######A) datos repartidos   
    #Segement your data by fold using the which() function 
    M_l_testing <- Base_mixed[i, ]
    M_l_training <- Base_mixed[-i, ]
    #Use the test and train data partitions however you desire...
    
    ##############
    ##############
    ##############Training
    ##############
    ###############
    #######A) para cada ps elijo que indices de ese tampnio estan mas cerca del PCA total para ambos
    #####en indices[1,] esta PCA MIX, en indices[2,] esta el PCA comun
    
    
    
    
    indice_final_PCA_glm[i,]=PCA_ambos(M_l_training,N0,.1,p_dim,ps_optimo_PCA_glm)[2,]
    
    
    indice_final_PCA_MIX_glm[i,]=PCA_ambos(M_l_training,N0,.1,p_dim,ps_optimo_PCA_mix_glm)[1,]
    
    
    #######B) pwe need to find the lienar combinations for that we need to accomoateh the variables
    
    
    
    #####PCA glm
    
    indices=indice_final_PCA_glm[i,]
    
    
    acomodate_data=accomodate_variables_for_PCA_PCAmix(indices,c,b,o,M_l_training)
  
    
    Total_final_PCA_glm=acomodate_data$XXsubset
    
    ######### #####index for PCAmix flm
    
    indices=indice_final_PCA_MIX_glm[i,]
    
    
    acomodate_data=accomodate_variables_for_PCA_PCAmix(indices,c,b,o,M_l_training)
    
    X_final=acomodate_data$Xsubset
    Cual_final=acomodate_data$Cualsubset
  
    
    #if the dimension of the conitnuos is 1 i need a matrix instead of a vector
    if(!is.null(X_final) & is.null(dim(X_final)[1])){X_final=as.matrix(X_final)}
  
    
    X_final_PCA_mix_glm=X_final
    Cual_final_PCA_mix_glm=Cual_final
    
    
    ##############
    ##############en M$score y M_PCA$score estan las combinaciones lineales
    
    
    M_glm=PCAmix(X.quanti = X_final_PCA_mix_glm, X.quali =Cual_final_PCA_mix_glm, ndim = 2, rename.level = TRUE,
                 weight.col.quanti = NULL, weight.col.quali = NULL, graph = FALSE)
    
    
    M_PCA_glm=prcomp(Total_final_PCA_glm,rank=1)
    
    
    
    
    ##########C) aplico glmfit y glm a los scorse (combinaciones lineales de PCA) en el training        
    
    #####miro los scores de cada una en el training
    
    M_glm_scores=M_glm$scores[,1]
    
    
    M_PCA_scores_glm=predict(M_PCA_glm,Total_final_PCA_glm)[,1]
    
    
    #######fiteo la logistica y fiteo la lineal para cada uno

    
    #####PCA_mix glm
    datos_PCA_mix_glm=data.frame(XX=M_glm_scores,Y=M_l_training[,1] )
    
    glmfit_PCA_mix_glm <- glm(Y~XX, data=datos_PCA_mix_glm, family = 'binomial'(link = "logit"))
    
    cutoff_optimo_pca_mix = cutoff.optimo.minMSE(glmfit_PCA_mix_glm,M_l_training[,1],M_glm_scores,trh)
    
    #####PCA_mix glm
    datos_PCA_glm=data.frame(XX=M_PCA_scores_glm,Y=M_l_training[,1] )
    
    glmfit_PCA_glm <- glm(Y~XX, data=datos_PCA_glm, family = 'binomial'(link = "logit"))
    
    cutoff_optimo_pca = cutoff.optimo.minMSE(glmfit_PCA_glm,M_l_training[,1],M_PCA_scores_glm,trh) 
    
    
    ##############
    ##############
    ##############TEST
    ##############
    ###############
    #######A)  le aplico esas combinaciones lineales al test para saber cuales son en el test
    
    #####acododo los datos en test porque sino el predict no adna
    
    
    #####acomodo los indices para PCA comun para GLM
    
    
    
    #######B) pwe need to find the lienar combinations for that we need to accomoateh the variables
    
    
    
    #####PCA glm
    
    indices=indice_final_PCA_glm[i,]
    
    
    
    acomodate_data=accomodate_variables_for_PCA_PCAmix(indices,c,b,o,M_l_testing)
    
    
    
    Total_final_PCA_glm_test=acomodate_data$XXsubset
    
    
    ######### #####index for PCAmix flm
    
    indices=indice_final_PCA_MIX_glm[i,]
    
    
    
    acomodate_data=accomodate_variables_for_PCA_PCAmix(indices,c,b,o,M_l_testing)
    
    X_final=acomodate_data$Xsubset
    Cual_final=acomodate_data$Cualsubset
    aux = names(Cual_final)
    
    #if the dimension of the conitnuos is 1 i need a matrix instead of a vector
    if(!is.null(X_final) & is.null(dim(X_final)[1])){X_final=as.matrix(X_final)}
    if(!is.null(Cual_final)){Cual_final=matrix(Cual_final, nrow = 1)
    colnames(Cual_final) <- aux}
    
    X_final_PCA_mix_glm_testing=X_final
    Cual_final_PCA_mix_glm_testing=Cual_final
    
    
    #############################
    #####predict  glm_MIX
    
    
    Mnew_glm=predict(M_glm,X.quanti= X_final_PCA_mix_glm_testing,X.quali=Cual_final_PCA_mix_glm_testing)
    
    #####predict PCA_glm
    
    
    aux=setdiff(names(Total_final_PCA_glm),names(Total_final_PCA_glm_test))
    if (length(aux)!=0){
      Total_final_PCA_glm_test[aux]=0}
    
    
    
    
    Mnew_PCA_glm=predict(M_PCA_glm, Total_final_PCA_glm_test)
    
    
    
    
    
    
    ##########################
    #######
    
    ##########B) predicgo con esso modelos en los scores aplicados al test
    # PCA mix
    Mnew_datos = data.frame(XX=Mnew_glm[,1],Y=M_l_testing[,1] )
    
    Mnew_aux_PCA_mix_glm = predict(glmfit_PCA_mix_glm, newdata =  Mnew_datos, type = "response") 
    
    yhat_PCA_mix = 1*(Mnew_aux_PCA_mix_glm > cutoff_optimo_pca_mix)
    
    costos_PCA_MIX[i] = 1-mean(yhat_PCA_mix == M_l_testing[,1])
    
    
    # PCA
    Mnew_datos = data.frame(XX=Mnew_PCA_glm[,1],Y=M_l_testing[,1] )
    
    Mnew_aux_PCA_glm = predict(glmfit_PCA_glm, newdata =  Mnew_datos, type = "response")
    
    yhat_PCA = 1*(Mnew_aux_PCA_glm > cutoff_optimo_pca)
    
    costos_PCA[i] = 1-mean(yhat_PCA == M_l_testing[,1])
    
    
    
    
    
    #########C) guardo esas predicciones para encontrear el corte optimo
    
    # resultados_PCA_MIX_glm = data.frame("fiteados"= Mnew_aux_PCA_mix_glm, "verdaderos"=M_l_testing[,1])
    # 
    # 
    # resultados_PCA_glm = data.frame("fiteados"= Mnew_aux_PCA_glm, "verdaderos"=M_l_testing[,1])
    # 
    
    
    
  
    
    
    
    
    #######D) miro para cada test i la funcion que estoy mirando que hace lo mejor posible
    
    # aux_pca_mix=cutoff.optimo(resultados_PCA_MIX_glm$verdaderos, resultados_PCA_MIX_glm$fiteados) 
    # #costo para discretas
    # costos_PCA_MIX[i]=aux_pca_mix[2]
    # 
    # 
    # aux_pca=cutoff.optimo(resultados_PCA_glm$verdaderos, resultados_PCA_glm$fiteados) 
    # 
    # costos_PCA[i]=aux_pca[2]
    
    
    
    #####confusion table for mix
    
    #corte_PCA_mix[i]=aux_pca_mix[1]
    
    
    # Mnew_compare = 1*(Mnew_aux_PCA_mix_glm>aux_pca_mix[1])
    # 
    # ct <- table(M_l_testing[,1], Mnew_compare)
    # 
    # 
   # ct <- table(M_l_testing[,1], yhat_PCA_mix)
    
    ######errror a reportar
    
   # tabla_pca_mix_glm[i,]=c(ct)
    
    
    
    #####confusion table for PCA
    
    # corte_PCA[i]=aux_pca[1]
    # 
    # Mnew_compare = 1*(Mnew_aux_PCA_glm>aux_pca[1])
    # 
    # ct <- table(M_l_testing[,1], Mnew_compare)
    
   # ct <- table(M_l_testing[,1], yhat_PCA)
    
    ######errror a reportar
    
   # tabla_pca_glm[i,]=c(ct)
    
    
    
    # if(i==1){salvar_glm_mix=cbind(rep(i,table(folds)[i]),Mnew_glm[,1],Mnew_aux_PCA_mix_glm,M_l_testing[,1])
    # salvar_glm_PCA=cbind(rep(i,table(folds)[i]),Mnew_PCA_glm[,1],Mnew_aux_PCA_glm,M_l_testing[,1])
    # 
    # }
    # 
    # if(i!=1){
    #   salvar_glm_mix=rbind(salvar_glm_mix,cbind(rep(i,table(folds)[i]),Mnew_glm[,1],Mnew_aux_PCA_mix_glm,M_l_testing[,1]))
    #   
    #   
    #   
    #   
    #   salvar_glm_PCA=rbind(salvar_glm_PCA,cbind(rep(i,table(folds)[i]),Mnew_PCA_glm[,1],Mnew_aux_PCA_glm,M_l_testing[,1]))
    #   
    # }
    
  }
  
  
  
  #####salvamos uno atras de otro los resultados de cross validation
  #####primera columna que i de cross validation (de 1 a 5)
  #####segunda columna los scores de pca para los testing
  #####tercera columna los fiteados de pca para testing
  #####cuarta columna los verdaderso testing
  
  # colnames(salvar_glm_mix)=c("scores","fiteados","verdaderos")
  # colnames(salvar_glm_PCA)=c("scores","fiteados","verdaderos")
  
  #write.table(salvar_glm_mix,'poverty_PCAmix_Pampeana_samequarter.txt',row.names=FALSE)
  #write.table(salvar_glm_PCA,'poverty_PCAmix_Pampeana_samequarter.txt',row.names=FALSE)
  
  
  
  #####printing the indeces
  
  freqindices=table(indice_final_PCA_MIX_glm)
  
  indice_final_glm_mix=as.numeric(names(freqindices)[order(freqindices, decreasing = TRUE)])[1:ps_optimo_PCA_mix_glm]
  
  
  
  
  freqindices=table(indice_final_PCA_glm)
  
  indice_final_glm_PCA=as.numeric(names(freqindices)[order(freqindices, decreasing = TRUE)])[1:ps_optimo_PCA_glm]
  
  # 
  # 
  # 
  # #####original variables
  # aux=c(c,b,o)
  v_o_e_glm_mix=intersect(aux,indice_final_glm_mix)
  #  
  v_o_e_glm_PCA=intersect(aux,indice_final_glm_PCA)
  # 
  # 
  # #######name of selected varibles
  # 
  # print("Choosen variables for PCA for poverty")
  # 
  # v_o_e_glm_PCA
  # 
  # 
  # 
  # print("Choosen variables for PCAmix for poverty")
  # 
  # v_o_e_glm_mix
  # 
  # 
  # ########tabla cross validataion para PCA_mix_glm
  # 
  # 
  # print("confussion table for PCAmix for poverty")
  # 
  # #tabla_pca_mix_glm
  # 
  # mat_PCA_mix = matrix(apply(tabla_pca_mix_glm, 2, mean), ncol = 2)
  # 
  # # TPR (sensitivity)
  # mat_PCA_mix[2,2]/(sum(mat_PCA_mix[2,]))
  # 
  # # TNR (Specificity)
  # mat_PCA_mix[1,1]/(sum(mat_PCA_mix[1,]))
  # 
  # # Misclass
  # (mat_PCA_mix[1,2]+mat_PCA_mix[2,1])/sum(mat_PCA_mix)
  # 
  # print("confussion table for PCA for poverty")
  # 
  # #tabla_pca_glm
  # mat_PCA = matrix(apply(tabla_pca_glm, 2, mean), ncol = 2)
  # 
  # # TPR (sensitivity)
  # mat_PCA[2,2]/(sum(mat_PCA[2,]))
  # 
  # # TNR (Specificity)
  # mat_PCA[1,1]/(sum(mat_PCA[1,]))
  # 
  # # misclass
  # (mat_PCA[1,2]+mat_PCA[2,1])/sum(mat_PCA)
  
  
  ###########
  
  indice_final_glm_mix=v_o_e_glm_mix
  
  indice_final_glm_PCA=v_o_e_glm_PCA
  
  MSE_PCA = mean(costos_PCA)
  MSE_PCA_mix = mean(costos_PCA_MIX)

  sd_PCA = sd(costos_PCA)
  sd_PCA_mix = sd(costos_PCA_MIX)
  
  return(list(MSE_PCA = MSE_PCA, sd_PCA = sd_PCA, MSE_PCA_mix = MSE_PCA_mix, sd_PCA_mix = sd_PCA_mix))
  
}


