function [Topt,Tsubopt] = getT4mix(X,H)
%Para el estadistico suficiente, no lo centramos (siguiendo el paper). Dejo
%comentado los centrados
[nh,q] = size(H);
[nx,p] = size(X);
Xc = centering(X);
Hc = centering(H);
if nx~=nh,
    error('number of rows in X must equate the number of rows in H');
end
k = q*(q+1)/2;
Hraro= zeros(nh,k);
for j=1:nh
    auxt = H(j,:)'*H(j,:);
    Hraro(j,:) = [diag(auxt);auxt(find(tril(ones(q,q),-1)))];
end
Topt = [X,Hraro];
%Topt = [Xc,Hraro];
Tsubopt = [X,H,Hraro];
%Tsubopt = [Xc,Hc,Hraro];
end