% Function to update the probability matrix when a gamma(r,j,k) is updated
%
% Call   : prob2 = update_prob(gnew,gold,j,k,r,x,y,prob)
%
% Arguments
% gnew   : new value of gamma(r,j,k)
% gold   : old value of gamma(r,j,k)
% j      : 2nd co-ordinate of the updated gamma parameter
% k      : 3rd co-ordinate of the updated gamma parameter
% r      : 1st co-ordinate of the updated gamma parameter
% x      : n x p covariate matrix
% y      : n x q response matrix
% prob   : n x q old probability matrix
%
% Values
% prob2  : n x q updated probability matrix 

function[prob2] = update_prob(gnew,gold,j,k,r,x,y,prob)

n = size(x,1);
if (r==0)
    xx = ones(n,1);
else
    xx = x(:,r);
end

prob2 = prob;

if (j ~= k)
deltaj  = (gnew-gold)* xx.*y(:,k);
deltak  = (gnew-gold)* xx.*y(:,j);


oldj = (1-prob(:,j))./prob(:,j);
oldk = (1-prob(:,k))./prob(:,k);

oldj(isnan(oldj)== 1) = 10^4;
oldk(isnan(oldk)== 1) = 10^4;

prob2(:,j) = 1./(1+oldj.*exp(-deltaj));
prob2(:,k) = 1./(1+oldk.*exp(-deltak));

else
    deltaj = (gnew-gold)*xx;
    oldj = (1-prob(:,j))./prob(:,j);
    oldj(isnan(oldj)== 1) = 10^4;
    
    prob2(:,j) = 1./(1+oldj.*exp(-deltaj));
end
end




