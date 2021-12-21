function JC=JqCqmatrix(q)

Aux =Cqmatrix(q);
Aux(any(Aux==1,2),:)=[];
JC=Aux;

