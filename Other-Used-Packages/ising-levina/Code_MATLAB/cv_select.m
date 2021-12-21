% function to choose the tuning parameter by Cross Validation for repeating
% experiments. 

function[lambda_opt] = cv_select(data_y, data_x, lambda, nfold, option)
K = size(data_y, 3);
lambda_opt = zeros(K,1);
matlabpool open;
%parpool open;
parfor k = 1:K
    sprintf('CV rep %d ', k)
    lambda_opt(k) = CVtune(data_y(:,:,k), data_x(:,:,k), lambda, nfold, option);
end
matlabpool close;
%parpool close;
end