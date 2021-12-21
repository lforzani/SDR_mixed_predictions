function [M,N] = PFC1(y,x)

[n,p] = size(x);
xo=x-ones(n,1)*mean(x);
N=1/n*(xo'*xo); 


yy=[y,y.^2];

yf=yy-ones(n,1)*mean(yy);
Pf = yf*inv(yf'*yf)*yf';
Xf = Pf*xo;
Sigmaf = 1/n*Xf'*Xf;

M =Sigmaf ;