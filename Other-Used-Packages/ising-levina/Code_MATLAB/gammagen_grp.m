% Function to generate parameter array for group lasso type of data  
%
% 
% Call    : gamma = gammagen_grp
%
% Arguments
% adj : q by q adjacent matrix
% p   : dimension of x
% val : the value of gamma(p,j,k)
% m   : number of effective predictors
% Values
% gamma: p+1 by q by q array with each gamma(:,j,k) being a vector with the
% fist m elements being non-zero with value val

function[gamma] = gammagen_grp(adj, p, val, m, prop)
q = size(adj,1);
gamma = zeros(p+1,q,q);
for j = 1:q
    for k = j:q
    if (k == j)
        gamma(2:p+1,j,j) = val* binornd(1, prop , p,1).*(2*binornd(1,0.5, p,1)-1);
        gamma(1,j,j) = val*(2*binornd(1,0.5)-1);
%         gamma(:,j,j) = val*binornd(1, prop, p+1,1);
%         gamma(1,j,j) = val;
    elseif (adj(j,k) == 1)
        gamma(2:(m+1),k,j) = val* (2*binornd(1,0.5,m,1)-1);
        gamma(1,k,j) = val* binornd(1,prop)*(2*binornd(1,0.5)-1);
        %gamma(:,k,j) = val*binornd(1, prop, p+1,1);
        gamma(:,j,k) = gamma(:,k,j);
    end
    end
end
end

        
        
        

