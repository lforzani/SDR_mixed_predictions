% function to calculate the L2-distance between a fitted parameter array
% and the true parameter.
% Arguments:
% gamma : true parameter array (p+1,q,q)
% gamma_mat: fitted parameter array (p+1, q, q, L)
% intercept: whether the L2 distance inlucde intercept. if 1 yes. if 0 no.
function[dist] = DistL2(gamma, gamma_mat, intercept)
p = size(gamma,1)-1;
q = size(gamma,2);
L = size(gamma_mat,4);
dist = zeros(L,1);
npar = (p+1)*q*(q+1)/2;

for j = 1:q
    for k = j:q
    tmp1 = gamma(:,k,j);
    tmp1 = tmp1(:,:,ones(L,1));
    tmp2 = reshape(gamma_mat(:,k,j,:), p+1, 1, L);
    diff = tmp1 - tmp2;
    if (intercept==0 && k==j)
        diff(1,:,:) = 0;
    end    
    dist = dist + reshape(sum(diff.^2,1), L,1);
    end
end
dist = dist./npar;
end