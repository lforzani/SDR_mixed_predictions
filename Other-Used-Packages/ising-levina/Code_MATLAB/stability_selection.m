%% Stability Selection for the Ising model with covariates
% INPUT:
% y : n*q response matrix
% x : n*p design matrix
% lambda: L*1 regularization parameters
% m: proportion of subsample as in the data
% K: number of replications
% Output:
% select_prob: the proportion of x being selected for each pairwise
% association and each lambda q*q*p*L
% select_max: the maximum selecting probability over lambda q*q*p


function[select_prob] = stability_selection(y,x,lambda, m, K,alpha, beta, omega)
[n q]= size(y);
p = size(x,2);
L = length(lambda);
n_sub = floor(n*m);
ind_gamma = zeros(p+1,q,q,L); %% record the selection of each 
                             % x_i over the repetitions for each lambda 
matlabpool open;
parfor k = 1:K
    k
    rng('default');
    rng(k);
    ind = randsample(n,n_sub);
    ytrain = y(ind,:);
    while(sum(sum(ytrain,1)==n_sub)>0 || sum(sum(ytrain,1)==0)>0)
        sprintf('resample')
    ind = randsample(n,n_sub);
    ytrain = y(ind,:);
    end
    gamma_temp = SolPathIdv(y, x, ind, lambda, alpha, beta, omega);
    ind_gamma = ind_gamma+(gamma_temp~=0);
    name = sprintf('DataAnalysis/New_stdX-4_stdY_nonstdLasso_gamma_mat%d.mat',k);
    iSaveGamma(name, gamma_temp);
    iClear(gamma_temp);
end
matlabpool close;

select_prob = zeros(q,q,p+1,L);
for l = 1:L
    for i = 1:p+1
    select_prob(:,:,i,l) = reshape(ind_gamma(i,:,:,l),q,q);
    end
end
select_prob = select_prob/K;
clear ind_gamma;
end


    
    

    













