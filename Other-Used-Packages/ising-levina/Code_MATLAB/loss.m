% Function to calculate the loss =  negative log-likelihood + l1 penalty
% 
% Call    : l = loss(y,x,par,lambda)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% par     : (p+1) x q x q parameter array
% lambda  : penalty parameter
%
% Values
% l       : scalar loss
function[loss neglikeli] = loss(y,x,par,lambda)

prob = cdprob(x,y,par);
ll = (prob.^y).*((1-prob).^(1-y));
l2 = log(ll);
s1 = sum(sum(l2));
if (isinf(s1)==1)
    s1 = -10^4;
end

q = size(y,2);
p = size(x,2);
n = size(y,1);

s2 = 0;
par1 = abs(par);

for m = 1:q
    v1 = sum(par1(2:p+1,m,m));
    v2 = sum(sum(par1(:,(m+1):q,m)));
    s2 = s2+ v1 + v2;
end

loss = -s1/(n)+lambda*s2;
neglikeli = -s1/n;
end
