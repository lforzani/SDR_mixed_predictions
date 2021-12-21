function [T,redu,proj,Ahat,Deltahat,fycent] = EM4mixture_continua(X,Y,g)
%UNTITLED Summary of this function goes here
%   INPUTS:
%  
%   X: continuas (nxp)
%  
%   Y: respuesta (nx1)
%   g: grado de la base para fy, si Y es continua, o 'disc'.
%
%   OUTPUTS:
%   redu: reduccion de tamanio r sin redudir
%   proj: datos proyectados

%   Modelo   X|Y = A fycent + error con varianza Delta

%   T = Xc
% ========================================================================
if ~ isstr(g)
    fycent = get_fy(Y,g);
else
    fycent = get_fyZ(Y);
end

n = size(X,1);
 
r = size(fycent,2);
 
% centra X
Xc = centering(X);

% analizamos X|Y
Pred = fycent;

% Coefficient of the regression for X|fy 
% OLS sin rango reducido Esto es A
Ahat = Xc'*Pred*inv(Pred'*Pred);

 
%estimador de . la covarianza de los residuos 

Deltahat = cov(Xc-Pred*Ahat');

% reducciones
redu = inv(Deltahat)*Ahat;
 

% --------estadistico suficiente
 
T = Xc;

% proyecciones
proj=T*redu;

 
