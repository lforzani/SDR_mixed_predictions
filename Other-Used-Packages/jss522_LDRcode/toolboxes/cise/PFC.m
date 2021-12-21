function [M,N] = PFC(y,x)

[n,p] = size(x);
xo=x-ones(n,1)*mean(x);
N=1/n*(xo'*xo); 


yy=[sqrt(abs(y)),y,y.^2];

yf=yy-ones(n,1)*mean(yy);
Pf = yf*inv(yf'*yf)*yf';
Xf = Pf*xo;
Sigmaf = (1/n)*Xf'*Xf;

M =(Sigmaf + Sigmaf')/2 ;
%M = Sigmaf;