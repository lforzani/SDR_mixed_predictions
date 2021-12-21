function [M,N] = SIR(y,x,nslice)

[n, p] = size(x);
xo=x-ones(n,1)*mean(x);
N=1/n*(xo'*xo); 


[a, index] = sort(y);
xmean = mean(x);
group = sort(mod(1:n, nslice)); 

sigmaeta = zeros(p, p);
diffmean = zeros(1, p);

for i=0:(nslice-1)
    diffmean = mean( x(index(group==i), :) ) - xmean;
    sigmaeta = sigmaeta + mean(group == i) * diffmean' * diffmean;
end

sigmaeta = (sigmaeta + sigmaeta') / 2;  

M = sigmaeta;
