function [Delta,A,beta] = get_initval_v2(Y,X,H,g)

[n,p] = size(X);
X = X - repmat(mean(X),n,1);

if ~isstr(g),
    fy = get_fy(Y,g);
else
    fy = get_fyZ(Y);
end

Pred = [fy H];
 

r=size(fy,2)

% Coefficient of the regression for X|H,fy
% OLS sin rango reducido
C = X'*Pred*inv(Pred'*Pred);
 
c2 = size(C,2);
 
%Coefficient of fy for the regression of X|H,fy
A = C(:,1:r);
 
%Coefficient of Hc for the regression of X|H,fy
beta = C(:,(r+1):c2);
 
Delta = cov(X-Pred*C');
%aux = inv(Fy'*Fy);
%Qf = eye(n)-Fy*aux*Fy';
%Delta = (X'*Qf*X)/n;
% alpha = inv(Delta)*X'*Fy*aux;
% nueva version alpha=AB
 %   ideltahalf = invsqrtm(Delta);
  %  [v,~] = firsteigs(ideltahalf*cov(X)*ideltahalf,dim);
   % A = ideltahalf*v;
    %B = inv(A'*Delta*A)*A'*X'*Fy*aux;
    %alpha = A*B;


