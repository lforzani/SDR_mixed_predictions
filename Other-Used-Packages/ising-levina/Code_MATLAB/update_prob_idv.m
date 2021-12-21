% Function to update the probability matrix when a gamma(r,j,k) is updated
%
% Call   : prob2 = update_prob(gnew,gold,j,k,r,x,y,prob)
%
% Arguments
% gnew   : new value of gamma(r,j,k)
% gold   : old value of gamma(r,j,k)
% j      : 2nd index of the updated gamma parameter
% k      : 3rd index of the updated gamma parameter
% r      : 1st index of the updated gamma parameter
% x      : n x p covariate matrix
% y      : n x q response matrix
% prob   : n x 1 old probability column
%
% Values
% prob_new  : n x 1 updated probability matrix 

function[prob_new] = update_prob_idv(gnew,gold,j,k,r,x,y, probj)

n = size(x,1);
    if (r==0)
    xx = ones(n,1);
    else
    xx = x(:,r);
    end

if (j ~= k)
    deltaj = (gnew-gold)*xx.*y(:,k);
else
    deltaj = (gnew-gold)*xx;
end
    oldj = (1-probj)./probj;
    oldj(isnan(oldj)== 1) = 10^4;
    prob_new = 1./(1+oldj.*exp(-deltaj));
end



