% Function to calculate the loss =  negative log-likelihood + l1 penalty
% 
% Call    : l = loss(y,x,par,lambda)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% par     : (p+1) x q parameter matrix
% lambda  : penalty parameter
%
% Values
% l       : scalar loss
function[loss likeli] = loss_idv(y,x,j,par,lambda)
n = size(y,1);
prob = cdprob(x,y,par);
ll = (prob(:,j).^y(:,j)).*((1-prob(:,j)).^(1-y(:,j)));
l2 = log(ll);
s1 = sum(l2);
if (isinf(s1)==1)
    s1 = -10^4;
end

s2 = sum(sum(abs(par(:,:,j))))-abs(par(1,j,j));

loss = -s1/(n)+lambda*s2;
likeli = -s1/n;
end
