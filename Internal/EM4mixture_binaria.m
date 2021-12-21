function [T, tau0,redu,proj,fycent] = EM4mixture_binaria(H,Y,g)
%UNTITLED Summary of this function goes here
%   INPUTS:
%    
%   H: binarias
%   Y: respuesta

% g =  grado de la base para fy, si Y es continua, o 'disc'.
%
%   OUTPUTS:
%   MODELO= H|Y los theta parametro naturales de H|Y son lineales 
%   en fycent y tienen coeficientes redu
%   
%   proj: datos proyectados
%   redu: matriz cuyas columnas son los vech de las matrices de
%   coeficientes del modelo Ising para cada r
%   
%   T: estadistico suficiente
%   T: estadistico suficiente ordenado segun Tauhatraro
%   T y redu estan acomodados primero los de la diagonal y luego el resto
% ========================================================================
 
       



if ~isstr(g),
    fycent = get_fy(Y,g);
else
    fycent = get_fyZ(Y);
end
  
[n,q] = size(H);
Hc = centering(H);
r=size(fycent,2);

 

% analizamos H|Y~BERNOULLI
tauhat = levina4ising(H,fycent,0);
tau0pre=tauhat(1,:,:);

k = q*(q+1)/2;
tau0=zeros(k,1);

auxx=reshape(tau0pre(1,:,:),[q,q]);
    auxx2= ones(q,q);
    tau0 = [diag(auxx);auxx(find(tril(auxx2,-1)))];

tauhat(1,:,:) = [];

  
 



  
%redu   tiene primero los correspondientes a la
%diagonal es decir a la parte de reduccion de H
%y luego la parte de la upper triangular que corresponde
%a la parte de arriba de la HH^T

redu= zeros(k,r);
for i = 1:r
    auxx=reshape(tauhat(i,:,:),[q,q]);
    auxx2= ones(q,q);
    redu(:,i) = [diag(auxx);auxx(find(tril(auxx2,-1)))];
end



 
  
% --------estadistico suficiente
HHt = zeros(n,q,q);
for j=1:n,
    HHt(j,:,:) = Hc(j,:)'*Hc(j,:);
end
 



%T, acomodo los H primero los H y luego los de
%upper triangular

T= zeros(n,k);
for j=1:n
    auxt = H(j,:)'*H(j,:);
T(j,:) = [diag(auxt);auxt(find(tril(ones(q,q),-1)))];
end


% proyecciones
proj=T*redu;
