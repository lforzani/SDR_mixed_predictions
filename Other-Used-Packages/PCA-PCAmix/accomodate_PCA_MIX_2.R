#######accomodate variables for PCA_mix

#c_aca =continuoas indeces
#b_aca=binary inidces
#o_aca=ordinal indices
#datos=data

###return : return(list(XXsubset = XXsubset,
######Xsubset=Xsubset continuas
#####Cual

accomodate_PCA_MIX_2<-function(c_aca,b_aca,o_aca,datos){


if (length(c_aca)==0){Xsubset=NULL}
if (length(c_aca)!=0){Xsubset=cbind(datos[,c_aca])}

if (length(c_aca)==0){X_final=NULL}
if (length(c_aca)!=0){X_final=cbind(datos[,c_aca])}


if (length(b_aca)==0){Hsubset=NULL}
if (length(b_aca)!=0){Hsubset=datos[,b_aca]}

if (length(o_aca)==0){Osubset=NULL}
if (length(o_aca)!=0){Osubset=datos[,o_aca]}



if(length(o_aca)==0 & length(b_aca)==0)
{Cual_final=NULL}


if(length(o_aca)==0 & length(b_aca)!=0) 
{Cual_final=apply(cbind(Hsubset), 2, factor)}




if(length(o_aca)!=0 & length(b_aca)==0)
{Cual_final=apply(cbind(Osubset), 2, factor)}



if(length(o_aca)!=0 & length(b_aca)!=0)
{Cual_final=apply(cbind(Hsubset,Osubset), 2, factor)}


  return(list(Xsubset=X_final , Cualsubset=Cual_final))
}
