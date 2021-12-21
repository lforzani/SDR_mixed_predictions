function [redu_optimal,redu_suboptimal, fycent,T_optimal,T_suboptimal,Tauhatraro0, Tauhatraro, Ahat, betahat, Deltahat] = EM4mixture_continua_binaria(X,H,Y,g)
%UNTITLED Summary of this function goes here
%   INPUTS:
%  
%   X: continuas (nxp)
%   H: binarias (nxq)
%   Y: respuesta (nx1)
%   g: grado de la base para fy, si Y es continua, o 'disc'.
%
%   OUTPUTS:
%   redu: base para span(eta-bareta) usando Levina
%   reduBL: reducion usando BIN laden
%   proj: datos proyectados
%   projBL: datos . proyectados Bin Laden
%   MODELO      X| H,  fycent   = A fycent + beta H
%       H| fycent = los parametros naturales son lineales
%                en fycent con Tauhat fy (Tauhat con levina
%       y Tauhatraro con Bin Laden
% ========================================================================
if ~ isstr(g)
    fycent = get_fy(Y,g);
else
    fycent = get_fyZ(Y);
end

n = size(X,1);
q = size(H,2);
r = size(fycent,2);
k = q*(q+1)/2;

Xc = centering(X);
Hc = centering(H);

% analizamos H|Y~BERNOULLI
tauhat = levina4ising(H,fycent);


tau0pre=tauhat(1,:,:);

k = q*(q+1)/2;
tau0 = zeros(k,1);

auxx = reshape(tau0pre(1,:,:),[q,q]);
auxx2 = ones(q,q);
tau0 = [diag(auxx);auxx(find(tril(auxx2,-1)))];

Tauhatraro0 = tau0;    
    
tauhat(1,:,:) = [];
% q arreglos y cada uno es de orden rxq 

% analizamos X|H,Y
if fycent==zeros(n,1),
    Pred = [Hc];
else
Pred = [fycent Hc];
end
% Coefficient of the regression for X|H,fy
% OLS sin rango reducido
%C = Xc'*Pred*inv(Pred'*Pred);
C=mvregress(Pred,Xc);
C=C';

c2 = size(C,2);

%Coefficient of fy for the regression of X|H,fy
Ahat = C(:,1:r);

%Coefficient of Hc for the regression of X|H,fy
betahat = C(:,(r+1):c2);

%Deltahat = cov(X-Pred*C');
Deltahat=cov(Xc-Pred*C');
                             
%aux1 = betahat'*inv(Deltahat)*Ahat;

Tauhatraro = zeros(k,r);
for i = 1:r
    auxx = reshape(tauhat(i,:,:),[q,q]);
    auxx2 = ones(q,q);
    Tauhatraro(:,i) = [diag(auxx);auxx(find(tril(auxx2,-1)))];
end

% --------estadistico suficiente
% chequear si va centrado

% % normalizamos para ver si se arreglan las cosas
% Deltahat = Deltahat/norm(Deltahat);
% Ahat = Ahat/norm(Ahat);
% Tauhatraro0 = Tauhatraro0/norm(Tauhatraro0);
% Tauhatraro = Tauhatraro/norm(Tauhatraro);
% betahat = betahat/norm(betahat);
%     

% reducciones
primero = inv(Deltahat)*Ahat;
if fycent==zeros(n,1),
    segundo = Tauhatraro(1:q,:);
else
segundo = Tauhatraro(1:q,:) - betahat'*primero;  %Tauhat-L;
end
tercero = Tauhatraro(q+1:k,:);
redu = [primero; segundo; tercero];

if fycent==zeros(n,1),
    seg = zeros(q,1);
else
    seg = -betahat'*primero;
end

ter=Tauhatraro(1:q,:);
cuar=Tauhatraro(q+1:k,:);
Tauhat=Tauhatraro;

Hraro = zeros(n,k);
for j=1:n
    auxt = H(j,:)'*H(j,:);
    Hraro(j,:) = [diag(auxt);auxt(find(tril(ones(q,q),-1)))];
end

T = [Xc,Hraro];

% proyecciones
proj=T*redu;


%Tauhatraro = csvread('redBinLaden.csv',1,0);
%Tau0BL=csvread('tau0.csv',1,0);
%segundo = Tauhatraro(1:q,:) - betahat'*primero;
%tercero = Tauhatraro(q+1:k,:);

%reduBL=[primero; segundo; tercero];

  
% proyecciones
%projBL=T*reduBL;
               
T1=[Xc,Hc];
T2=Hraro;
seg=-betahat'*primero;
T_optimal=[Xc,Hraro];
T_suboptimal=[Xc,Hc,Hraro];

redu_optimal=[redu];
if fycent==zeros(n,1)
    redu_suboptimal =[primero; ter;cuar];
else
redu_suboptimal=[primero;seg;ter;cuar];
end
end



