
function[tmp] = get_gamma(beta, beta0, p, q, j)
tmp = reshape([beta0;beta], p+1,q);
tmp1 = tmp;
tmp(:,j) = tmp1(:,1);
tmp(:,1:j-1) = tmp1(:,2:j);
end