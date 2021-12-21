% Function to compute the entire solution path by decreasing lambda
% gradually starting from the maximum effective shrinkage parameter
%
% Call         : [par_path lambda_path a] = SolPath(y,x,lambda_max, lambda_min, k)
%
% Arguments
% y            : n x q response matrix
% x            : n x p covariate matrix
% lambda_max   : a vector of exponentially decreasing lambda
%
% Values
% par_path     : (p+1) x q x q x L array where each (:,:,:,i) the entry 
%                denotes the (p+1) x q x q parameter estimate at lambda(i)  


function[par_path] = SolPath(y,x,lambda,option)
tic;
p = size(x,2);
q = size(y,2);
L = length(lambda);
par_old = initialize(y,x);
%par_old = 2*ones(p+1,q,q);

par_path = zeros(p+1, q, q, L);
if strcmp(option, 'joint')
    for i = 1:length(lambda)
  %  sprintf('the %dth lambda', i)
    par_path(:,:,:,i) = CoordinateShooting(y,x,lambda(i),par_old);
    par_old = par_path(:,:,:,i);
    end
elseif strcmp(option, 'separate')
    for i = 1:length(lambda)
    %sprintf('the %dth lambda', i)
    par_path(:,:,:,i) = CoordinateShooting_idv(y,x,lambda(i),par_old);
    par_old = par_path(:,:,:,i);
    end
end
end
