% function to choose the tuning parameter by validating on another
% independent data set
% Arguments
% data_y: pool of y data (n,q,K)
% data_x: pool of x data (n,q,K)
% gamma_mat: pool of fitted parameters for each data and solution path...
%             (p+1, q, q, L, K)
function[lambda_opt likeli] = Validtune(data_y, data_x, gamma_mat)
L = size(gamma_mat, 4);
K = size(data_y, 3);
lambda_opt = zeros(K,1);
likeli = zeros(L,K);
for k = 1:K
    if (k==K)
        testInd = 1;
    else 
        testInd = k+1;
    end
    for l = 1:L
        likeli(l,k) = loss(data_y(:,:,testInd), data_x(:,:,testInd), gamma_mat(:,:,:,l,k), 0);
    end
    [c, lambda_opt(k)] = min(likeli(:,k));
    clear c;    
end
end

