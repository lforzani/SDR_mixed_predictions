% function to choose the best sensitivity and specificity for repeating
% experiments. 

function[select_prob] = stability_sim_select(data_y, data_x, lambda, option, reps)

[n q K]= size(data_y);
p = size(data_x,2);
L = length(lambda);
n_sub = n/2;
select_prob = zeros(p+1,q,q,L,K); %% record the selection freq of parameters.
matlabpool open;
parfor k = 1:K
for r = 1:reps
    sprintf('stability select reps%d round%d', k, r)
    ind = randsample(n,n_sub);
    ytrain = data_y(ind,:,k);
%     while(sum(sum(ytrain,1)==n_sub)>0 || sum(sum(ytrain,1)==0)>0)
%         sprintf('resample')
%     ind = randsample(n,n_sub);
%     ytrain = data_y(ind,:,k);
%     end
    xtrain = data_x(ind,:,k);
    gamma_temp = SolPath(ytrain, xtrain, lambda, option);
    select_prob(:,:,:,:,k) = select_prob(:,:,:,:,k)+(gamma_temp~=0);
end
end
matlabpool close;
select_prob = select_prob/reps;
select_prob = max(select_prob, [], 4);
select_prob = reshape(select_prob, p+1, q,q,K);
end

