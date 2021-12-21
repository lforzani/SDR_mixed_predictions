function [Y, H TrueRedBin] = GenDataBinconditional2(n, q, r)

%Generaci?n de datos para el caso de predictores binarios. Modelo para la regreci?n inversa (ising model) con
%H|Y = Gamma1 fy[1]+Gamma2 fy[2]

 
%q = 4; % cantidad de discretas
  
%r = 2; %dos predictores  
 
%n = 5000; % tamanio de muestra
 
%H|Y = Gamma1 fy[1]+Gamma2 fy[2]
% generamos coeficientes de la regresion de H en y
%Gamma es p por p simetrica porque el eta_12=eta_21 

 Gamma1 = zeros(q,q);
% Gamma1(1,1) = 1;
% Gamma1(2,2) = 1;
% Gamma1(3,3) = 1;
% Gamma1(4,4) = 1;
% Gamma1(1,2) = 10;
% Gamma1(2,1) = 10;
% Gamma1(1,3) = 0;
% Gamma1(3,1) = 0;
% Gamma1(3,4) = 10;
% Gamma1(4,3) = 10;
% Gamma1 = Gamma1/sqrt(204);
% Gamma1 = Gamma1*2;
% 
% Gamma2 = 2*Gamma1;


Gamma1(1,1) = 1;
Gamma1(2,2) = 1;
Gamma1(3,3) = 1;
Gamma1(4,4) = 1;
Gamma1(1,2) = 10;
Gamma1(2,1) = 10;
Gamma1(1,3) = 0;
Gamma1(3,1) = 0;
Gamma1(3,4) = 10;
Gamma1(4,3) = 10;


Gamma1 = Gamma1/sqrt(204);
Gamma1 = Gamma1*2;
Gamma1 = Gamma1*2;

Gamma2 = 2*Gamma1;
 
% Generaramos la respuesta y fY
g = r+1;
probs = [0.32, 0.28, 0.40];

Y =  randsample(g, n,true, probs);

fy = zeros(n,r);

xx = unique(Y); 
nk = zeros(size(xx)); 
for i = 1:length(xx) 
    nk(i) = sum(Y == xx(i)); 
end
    
   
for i = 1:n
  for j = 1:r
    
    if Y(i)==j
        IND = 1;
    else
        IND = 0;
    end
    
    fy(i,j) = IND  - nk(j)/n;
  end
end
 
 
  %CENTRAR fy
  meanfy = mean(fy, 1);
  fycent(:,1) = fy(:,1) - meanfy(:,1);
  fycent(:,2) = fy(:,2) - meanfy(:,2);
  
  % generamos las discretas
  
 %for i = 1:n 
  %  aux(i,:,:) = Gamma1*fycent(i,1)+Gamma2*fycent(i,2);
 % end
  
   
      auu1=unique(fycent,'rows');
  for s=1:(r+1)
   nn(s)=sum(ismember(fycent,auu1(s,:),'rows'));
   Ind(:,s)=ismember(fycent,auu1(s,:),'rows');
  end
  
  
  
  
    
   nreps=10^6;
  H = binornd(1,1/2,n,q); 
  
  
  
  for j=1:nreps;
  
  
   for s=1:(r+1);
        aux(s,:,:) = Gamma1*auu1(s,1)+Gamma2*auu1(s,2);
      aux1=reshape(aux(i,:,:),[q, q]);
      etass=aux1;
       
       
      
      for hh=1:q;
          
         Yraro=H(ismember(fycent,auu1(s,:),'rows'),:);
     
          
          Yraro(:,hh)=1;
          
          etas=Yraro*etass;
          
          ss=etas(hh);
          
          H(ismember(fycent,auu1(s,:),'rows'),hh)=binornd(1,1./(1+exp(-ss)),nn(s),1);
      end
  end
  
  end
  
  tauhat = levina4ising(H,fycent,0)
tau0pre=tauhat(1,:,:);
  
  
  
  auxx2 = ones(q,q);
  TrueRedBin = [diag(Gamma1); Gamma1(find(tril(auxx2,-1)))];
   
  