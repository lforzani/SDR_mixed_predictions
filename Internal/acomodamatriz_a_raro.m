function [T, redu,proj,Tauhat,fycent] = EM4mixture_binaria(H,Y,g)
%UNTITLED Summary of this function goes here
%   INPUTS:
%    
%   H: binarias
%   Y: respuesta

% g =  grado de la base para fy, si Y es continua, o 'disc'.
%
%   OUTPUTS:
%   redu: base para span(eta-bareta)
%   proj: datos proyectados
% ========================================================================
 
  


if ~isstr(g),
    fycent = get_fy(Y,g);
else
    fycent = get_fyZ(Y);
end
  
[n,q] = size(H);
Hc = centering(H);
r=size(fycent,2)

 

% analizamos H|Y~BERNOULLI
tauhat = levina4ising(H,fycent,0);
tauhat(1,:,:) = [];

  
 

k = q*(q+1)/2;

Tauhat = zeros(k,r);
L = zeros(k,r);
for i = 1:r
    Tauhat(:,i) = vech(reshape(tauhat(i,:,:),[q,q]));
end

% reducciones
 
 
redu = Tauhat;

  
% --------estadistico suficiente
HHt = zeros(n,q,q);
for j=1:n,
    HHt(j,:,:) = Hc(j,:)'*Hc(j,:);
end
T = vechmd2(HHt)';

% proyecciones
proj=T*redu;