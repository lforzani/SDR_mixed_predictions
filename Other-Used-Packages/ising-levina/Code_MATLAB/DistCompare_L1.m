
function[dist] = DistCompare_L1(gamma, gamma_mat)
% Input
% gamma: true parameter array, p+1*q*q
% gamma_mat: fitted path parameter aray, p+1*q*q
%
% Output
% dist : vector of L1 distance between gamma and member of gamma_mat

q = size(gamma,2);
L = size(gamma_mat,4);

dist = zeros(L,1);

for l = 1:L
    delta = abs(gamma-gamma_mat(:,:,:,l));
    for j = 1:q
        dist(l) = dist(l) + sum(sum(delta(:,j:q,j)));
    end
end

        
            
        

