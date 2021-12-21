% Function to calculate the loss =  negative log-likelihood + l1 + group
% penalty
% 
% Call    : l = loss(y,x,par,lambda1, lambda2)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% par     : (p+1) x q x q parameter array
% lambda1 : L1 penalty parameter
% lambda2 : Group-lasso penalty parameter
% Values
% L : loss
function[l] = loss_grp(y,x,par,lambda1, lambda2)

prob = cdprob(x,y,par);
l2 = log((prob.^y).*((1-prob).^(1-y)));
s1 = sum(sum(l2));

q = size(y,2);
p = size(x,2);
n = size(y,1);

s2 = 0;
par1 = abs(par);
square = zeros(p,1);

for m = 1:q
    v1 = sum(par1(2:p+1,m,m));
    v2 = sum(sum(par1(1,(m+1):q,m)));
    s2 = s2+ v1 + v2;
    square = square + sum(par1(2:p+1,m+1:q,m).^2 ,2);    
end
square = sqrt(square);
l = -s1/(n)+lambda1*s2 + lambda2/2*sum(square);
end
