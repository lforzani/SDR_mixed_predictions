###accomodate data for applying PCA and PCAmix

#inputs
###c are the columns that have continuous variables
###o are the columns that have ordinal variables
###b are the columns that have binary variables
###sub is a subset of 1:length(c(c,b,o))
###datos, matrix of data

#output
###Xsubset continuos part for PCAmix
###Cualsubset binary + ordinal part for PCAmix
###XXsubset (all for PCA)

accomodate_variables_for_PCA_PCAmix<-function(sub,c,b,o=NULL,datos){
  
  library(fastDummies)
  
  Xsubset=NULL
  Hsubset=NULL
  Osubset=NULL

  if (length(intersect(sub, c))==0){Xsubset=NULL}
  if (length(intersect(sub, c))!=0){Xsubset=datos[,c(intersect(sub, c))]}
  
  if (length(intersect(sub, b))==0){Hsubset=NULL}
  if (length(intersect(sub, b))!=0){Hsubset=datos[,c(intersect(sub, b))]}
  
  if (length(intersect(sub, o))==0){Osubset=NULL}
  if (length(intersect(sub, o))!=0){Osubset=datos[,c(intersect(sub, o))]}
  
 
    
  if(length(intersect(sub, o))==0 & length(intersect(sub, b))==0 & length(intersect(sub, c))!=0)
  {Cualsubset=NULL;
  XXsubset=cbind(Xsubset)}
  
  
  if(length(intersect(sub, o))==0 & length(intersect(sub, b))!=0 & length(intersect(sub, c))==0)
  {Cualsubset=apply(cbind(Hsubset), 2, factor);
  XXsubset=cbind(Hsubset)}
  
  
  if(length(intersect(sub, o))==0 & length(intersect(sub, b))!=0 & length(intersect(sub, c))!=0)
  {Cualsubset=apply(cbind(Hsubset), 2, factor);
  XXsubset=cbind(Xsubset,Hsubset)}
  
  
 
  if(!is.null(o)){
    if(length(intersect(sub, o))!=0 & length(intersect(sub, b))==0 & length(intersect(sub, c))==0)
    {Cualsubset=apply(cbind(Osubset), 2, factor);
    Osubset=cbind(Osubset)
    Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
    XXsubset=cbind(Osubset)}
    
    
    if(length(intersect(sub, o))!=0 & length(intersect(sub, b))==0  & length(intersect(sub, c))!=0)
    {Cualsubset=apply(cbind(Osubset), 2, factor);
    Osubset=cbind(Osubset)
    Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
    XXsubset=cbind(Xsubset,Osubset)}
    
    
    if(length(intersect(sub, o))!=0 & length(intersect(sub, b))!=0 & length(intersect(sub, c))!=0)
    {Cualsubset=apply(cbind(Hsubset,Osubset), 2, factor);
    Osubset=cbind(Osubset)
    Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
    XXsubset=cbind(Xsubset,Hsubset,Osubset_dummy)}
    
    
    if(length(intersect(sub, o))!=0 & length(intersect(sub, b))!=0 & length(intersect(sub, c))==0)
    {Cualsubset=apply(cbind(Hsubset,Osubset), 2, factor);
    Osubset=cbind(Osubset)
    Osubset_dummy=dummy_cols(Osubset,select_columns=c(names(Osubset)))[,-(1:dim(Osubset)[2])];
    XXsubset=cbind(Hsubset,Osubset_dummy)}
  
    }
  
  return(list(XXsubset = XXsubset,
              Xsubset=Xsubset , Cualsubset=Cualsubset))
  
}
