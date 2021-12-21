function T = getTbin(H)
[n,q] = size(H);
k = q*(q+1)/2;
T= zeros(n,k);
for j=1:n
    auxt = H(j,:)'*H(j,:);
T(j,:) = [diag(auxt);auxt(find(tril(ones(q,q),-1)))];
end
end