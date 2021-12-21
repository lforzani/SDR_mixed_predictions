% Function to initialize beta using penalized logistic regression
%
% Call    : a = initialbeta(y, x, lambda)
%
% Arguments
% y       : n x q response matrix
% x       : n x p covariate matrix
% lambda  : penalty parameter
%
% Values
% a       : returned scalar initial value

function [ a ] = initialbeta(y, x, lambda)

%% y is n_dim 0-1 vector, x is n_dim cts vector.
n = size(y,1);
diff = 1;
delta = 10^-6;
thres = 10^-3;

%% initialize
beta = [1 ; 1];
beta_new = beta;

p = cdprob(x,y,beta);
w = p.*(1-p);

while(diff>=delta)
    t1 = w'*(x.^2);          %% update coefficients
    t2 = x'*(y-p);
    if (abs(t1)>=thres)
    beta_new(2) = lasso(t1*beta(2)+t2, n*lambda)/t1;
    end
    
    if (abs(sum(w))>=thres)   %% update intercept
    beta_new(1) = sum(y-p)/sum(w)+beta(1);
    end
    
    diff = sum(abs(beta-beta_new));
    beta = beta_new;
    
    p = cdprob(x,y,beta);     %% update parameter functions
    w = p.*(1-p);
end

a = beta(2);

end

