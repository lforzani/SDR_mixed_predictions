function LC = LqCqmatrix(q)

Aux = Cqmatrix(q);
Aux(any(Aux==1/2,2),:)=[];
LC = Aux;

