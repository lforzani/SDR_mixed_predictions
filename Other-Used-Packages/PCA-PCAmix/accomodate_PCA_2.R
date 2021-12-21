#######accomodate variables for PCA

#c_aca =continuoas indeces
#b_aca=binary inidces
#o_aca=ordinal indices
#datos=data

###return : TOTAL FINAL PCA

accomodate_PCA_2<-function(c_aca,b_aca,o_aca,datos){
  
  if (length(c_aca)==0){Xsubset=NULL}
  if (length(c_aca)!=0){Xsubset=datos[,c_aca]}
  
  if (length(b_aca)==0){Hsubset=NULL}
  if (length(b_aca)!=0){Hsubset=datos[,b_aca]}
  
  if (length(o_aca)==0){Osubset=NULL}
  if (length(o_aca)!=0){Osubset=datos[,o_aca]}
  
  
  if(length(o_aca)==0 & length(b_aca)==0 & length(c_aca)!=0)
  {Total_final_PCA=cbind(Xsubset)}
  
  
  if(length(o_aca)==0 & length(b_aca)!=0 & length(c_aca)==0)
  {Total_final_PCA=cbind(Hsubset)}
  
  
  if(length(o_aca)==0 & length(b_aca)!=0 & length(c_aca)!=0)
  {Total_final_PCA=cbind(Xsubset,Hsubset)}
  
  if(length(o_aca)!=0 & length(b_aca)==0 & length(c_aca)==0)
  {Osubset=cbind(Osubset);
  Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
  Total_final_PCA=cbind(Osubset)}
  
  if(length(o_aca)!=0 & length(b_aca)==0  & length(c_aca)!=0)
  {Osubset=cbind(Osubset);
  Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
  Total_final_PCA=cbind(Xsubset,Osubset)}
  
  
  if(length(o_aca)!=0 & length(b_aca)!=0 & length(c_aca)!=0)
  {Osubset=cbind(Osubset);
  Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
  Total_final_PCA=cbind(Xsubset,Hsubset,Osubset_dummy)}
  
  
  if(length(o_aca)!=0 & length(b_aca)!=0 & length(c_aca)==0)
  {Osubset=cbind(Osubset);
  Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
  Total_final_PCA=cbind(Hsubset,Osubset_dummy)}
 
  
  
  Total_final_PCA 
}