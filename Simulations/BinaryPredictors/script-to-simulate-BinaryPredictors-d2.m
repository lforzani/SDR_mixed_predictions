clear all;
%Simulation for binaries predictor with d=2
n =2000;
q = 10; 
nreps=100;

k=[100 200 300 500 750];
r = 5;
d=2;

kk = q*(q+1)/2;

N = length(k);
nmax = max(k);

Gamma1 = zeros(q,q);
Gamma1(1,1) = 1;
Gamma1(2,2) = 1;
Gamma1(3,3) = 1;
Gamma1(4,4) = 1;
Gamma1(1,2) = 30;
Gamma1(2,1) = 30;
Gamma1(1,3) = 5;
Gamma1(3,1) = 5;
Gamma1(3,4) = 30;
Gamma1(4,3) = 30;
Gamma1(2,3) = 10;
Gamma1(3,2) = 10;
for i=2:q
    Gamma1(i,i)=1;
    Gamma1(i-1,i)=30;
    Gamma1(i,i-1)=30;
end


for i=7:q
Gamma1(i,i) = 0;
Gamma1(i-1,i) = 0;
Gamma1(i,i-1) = 0;
end

Gamma1=Gamma1/sqrt(sum(sum(Gamma1)));


Gamma2=zeros(q,q);


for i=1:6
Gamma2(i,i)=1;
end

Gamma2=Gamma2/sqrt(sum(sum(Gamma2)));

Gamma=3*[Gamma1,4*Gamma2,Gamma1,Gamma1,Gamma1];

Gannas=reshape(Gamma,[q,q,r]);

auxx2 = ones(q,q);
TrueRedBin=[diag(Gamma1) diag(Gamma2);Gamma1(find(tril(auxx2,-1))) Gamma2(find(tril(auxx2,-1)))];

TRB1 = [diag(Gamma1) ; Gamma1(find(tril(auxx2,-1)))];
TRB2 = [diag(Gamma2) ; Gamma2(find(tril(auxx2,-1)))];

%To count true selected variables
%true zeros
 TZ1 = find(TRB1==0);
 TZ2 = find(TRB2==0);
 %selected variables (non zeros)
 TNZ1 = find(TRB1~=0); 
 TNZ2 =find(TRB1~=0);
 
indexZ=zeros(q,q);
for i=1:q
    m=i;
    indexZ(i,1)=i;
    for j=2:i
    indexZ(i,j)=m+(q+1-j);
    m= indexZ(i,j);
    end 
    if(i<q)
    indexZ(i,i+1)=indexZ(i,i)+ (q-i+1);
    end
    for j=(i+2):q;
        aux =indexZ(i,j-1);
        indexZ(i,j) = aux +1;
    end
end

CoincidZ=zeros(nreps,length(k));
CoincidNZ =zeros(nreps, length(k));

TotalCoincidZ = zeros(nreps,length(k));
TotalCoincidNZ =zeros(nreps, length(k));

%%
g = r+1;
probs = ones(g,1)/g;
YY1 =  randsample(g, n,true, probs);
fycent1=get_fyZ(YY1);

[HH1 ] = GenDataBinariasOPv2(n, fycent1, Gannas);

Top=getTbin(HH1);
  
Proy_True = Top*TrueRedBin;

for jj =1:nreps;
    jj
    YY =  randsample(g, k(N),true, probs);
    
     fycent=get_fyZ(YY);
   
   
    [HH] = GenDataBinariasOPv2(k(N), fycent, Gannas);
    

    for mm = 1:length(k)
       Y=YY(1:k(mm));
       H=HH(1:k(mm),:);
       TT_optimal =Top(1:k(mm),:);
        
 
 %  Parameter estimation without reduction                
  [T, tau0,redu,proj,fycent] = EM4mixture_binaria(H,Y,'disc');
         
 % Optimal reduction with automatic variable selection
  [red_est] = rr4bin(H,Y,d,'disc','auto');
    if size(red_est,2) < d,
        red_est(:,d)=zeros(size(red_est,1),1);
    end
    
    Proy_est = Top*red_est;

red_est1 = red_est(:,1);
red_est2 =red_est(:,2);

for i=1:q
if (red_est1(indexZ(i,:))==0 & red_est2(indexZ(i,:))==0)
CountZeros(i) =1;
else
CountZeros(i)=0;
end
end
TrueZeros=sum(CountZeros(7:q)==1)/4;
FalseZeros=(sum(CountZeros(1:6)==1))/6;

CoincidZ(jj,mm) =TrueZeros;
CoincidNZ(jj,mm) =1-FalseZeros;

%To count coinciden in total zeros (no seleceted variables)
  TotalTrueZeros=(sum(red_est1(TZ1)==0) + sum(red_est2(TZ2)==0))/(length(red_est1(TZ1))+length(red_est2(TZ2)));
 TotalFalseZeros=1-((sum(red_est1(TNZ1)==0) + sum(red_est2(TNZ2)==0))/(length(red_est1(TNZ1))+length(red_est2(TNZ2))));

  
TotalCoincidZ(jj,mm) = TotalTrueZeros;
TotalCoincidNZ(jj,mm) =1-TotalFalseZeros;


%PCA Filmer & Pritchet over H and interactions   
[coefPCA, scoresPCA] = pca(TT_optimal);
reduPCA=coefPCA(:,1:d);
Proy_PCA = Top*reduPCA; 
       
%PFC taking all variables as continuous  
[Proy_est_PFCTT,redu_est_PFC] = ldr(Y,TT_optimal,'pfc','disc',d);
Proy_est_PFC = Top*redu_est_PFC;

% Compute distante (Frobenious and Eucliean) between true and estimated reduction 
[Roptimalauto_F(jj,mm) Roptimalauto_2(jj,mm)] = diferencia(red_est,TrueRedBin);

[R_PCA_F(jj,mm) R_PCA_2(jj,mm)] = diferencia(reduPCA,TrueRedBin);

[R_PFC_F(jj,mm) R_PFC_2(jj,mm)] = diferencia(redu_est_PFC,TrueRedBin);

% Compute distante (Frobenious and Eucliean) between out-of-sample projections

[PDoptimalauto_F(jj,mm) PDoptimalauto_2(jj,mm)]= diferencia(Proy_est,Proy_True);

[PD_PCA_F(jj,mm) PD_PCA_2(jj,mm)]= diferencia(Proy_PCA,Proy_True);

[PD_PFC_F(jj,mm) PD_PFC_2(jj,mm)]= diferencia(Proy_est_PFC,Proy_True);

    end
end

%%


cd '~/Figures/BinaryPredictors'

csvwrite('Binary-d2-estimator.csv',Roptimalauto_2);
csvwrite('Binary-d2-proyection.csv',PDoptimalauto_2);
          
csvwrite('Binary-d2-estimator_PCA.csv',RD_PCA_2);
csvwrite('Binary-d2-proyection_PCA.csv',PD_PCA_2);
 
csvwrite('Binary-d2-estimator_PFC.csv',RD_PFC_2);
csvwrite('Binary-d2-proyection_PFC.csv',PD_PFC_2);
